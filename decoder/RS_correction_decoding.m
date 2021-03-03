% Copyright (c) 2021 by Jae-Won Kim, Jaeho Jeong, and Seong-Joon Park,
% from Coding and Cryptography Lab (CCL), Department of Electrical and Computer Engineering,
% Seoul National University, South Korea.
% Email address: jaehoj@ccl.snu.ac.kr
% All Rights Reserved.

% Code for proposed decoding method with C4
% input files: science_raw_read%d.txt, science_raw_cnt%d.txt
% output file: decoding_final_result.txt
% image_restart.txt: used only for comparing the result

clear;

%%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LT_L = 256; %number of bits in payload
LT_seed_bit = 4*8;
RS_parity_num = 2;
LT_K = 16050;
LT_c = 0.025;
LT_delta = 0.001;
LT_s = LT_c * sqrt(LT_K) * log(LT_K/LT_delta);

RNG_input_seed = 10;

soliton = zeros(1,LT_K);
for i=1:LT_K
    if(i==1)
        soliton(i) = 1/LT_K;
    end
    
    if(i~=1)
        soliton(i) = 1/i/(i-1);
    end
end

tau = zeros(1,LT_K);
KSratio = round(LT_K / LT_s);
for i=1:LT_K
    if(i == KSratio)
        tau(i) = LT_s * log(LT_s/LT_delta) / LT_K;
    end
    
    if(i<KSratio)
        tau(i) = LT_s / LT_K / i;
    end
end

LT_Z = 0;
for i=1:LT_K
    LT_Z = LT_Z + soliton(i) + tau(i);
end

Robust_soliton = zeros(1,LT_K);
for i=1:LT_K
    Robust_soliton(i) = (soliton(i) + tau(i)) / LT_Z;
end

sum_test = zeros(1,LT_K);

for i=1:LT_K
    if(i == 1)
        sum_test(i) = Robust_soliton(1);
    else
        sum_test(i) = sum_test(i-1) + Robust_soliton(i);
    end
end

rng(RNG_input_seed);
Binary_Data_input = zeros(LT_K,LT_L);

input_filename = sprintf('image_restart.txt');
[FP] = fopen(input_filename,'r');

input_data_save = fscanf(FP,'%d');
fclose(FP);

input_file_bit_size = size(input_data_save,1);


for i=1:LT_K
    for j=1:LT_L
        if((i-1)*LT_L+j<input_file_bit_size + 1)
            Binary_Data_input(i,j)=input_data_save((i-1)*LT_L+j);
        else
            Binary_Data_input(i,j)=randi(2)-1;
        end     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%For modification%%%%%%%%%%%%%%%%%%%%%%%%
start_index = 1;
end_index = 6;
Result_collect = zeros(1,end_index-start_index+1);  %% 0 - lack of inference, 1 - decoding success, 2 - decoding fail
Trial_Number_collect = zeros(1,end_index-start_index+1);  %% marking the number of trials
LT_real_N_collect = zeros(1,end_index-start_index+1);
LT_decoding_observance_collect = zeros(1,end_index-start_index+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%LT parameter + seed parameter%%%%%%%%%%%%%%%%%%%
for recursive=start_index:end_index
    
    raw_input_filename = sprintf('raw_read%d.txt',recursive);
    raw_input_number = sprintf('raw_cnt%d.txt',recursive);          %% for cluster size 2
    [FP1] = fopen(raw_input_filename,'r');
    [FP2] = fopen(raw_input_number,'r');
    
    raw_data_save = fscanf(FP1,'%d');
    raw_number_save = fscanf(FP2,'%d');
    
    fclose(FP1);
    fclose(FP2);

    number_of_strands = size(raw_data_save,1)/(LT_seed_bit+LT_L+8*RS_parity_num);
    LT_parity_candidate = zeros(number_of_strands,LT_seed_bit+LT_L+8*RS_parity_num);

    for i=1:number_of_strands
        for j=1:(LT_seed_bit+LT_L+8*RS_parity_num)
            LT_parity_candidate(i,j)=raw_data_save((i-1)*(LT_seed_bit+LT_L+8*RS_parity_num)+j);
        end
    end


    LT_real_N = 0;         %during decoding, number of clusters saved for future inference (when degree one comes)
    LT_decoding_observance = 0; %number of strands used for the decoding (count when RS ok but not inferred)


    Inferred_bitstream = zeros(LT_K,LT_L);
    Inferred_idx = zeros(LT_K,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    temp_vector_candidate = zeros(1,LT_seed_bit);
    
    corrected_strand = [];
    trial_number = 0;  % line number in input file
    real_trial = 0;    % trial number
    k = 0;
    while((real_trial < number_of_strands) && (sum(Inferred_idx)~=LT_K))
       flag=0;
        if (trial_number < number_of_strands)
            trial_number = trial_number + 1;
            
            RS_check_indicator = 0;
            RS_test_temp = zeros(1,(LT_seed_bit+LT_L+8*RS_parity_num)/8);
            eightbit_temp = zeros(1,8);
            binary_candidate = zeros(1,LT_seed_bit+LT_L+8*RS_parity_num);            
%%            
            if(raw_number_save(trial_number,1)==1)
                
                for i=1:(LT_seed_bit+LT_L+8*RS_parity_num)/8
                    for j=1:8
                        eightbit_temp(j)=LT_parity_candidate(trial_number,8*i-j+1);
                    end
                    RS_test_temp(i) = bi2de(eightbit_temp); 
                end
                
                RS_codeword_test = gf(RS_test_temp,8);
                [~,cnumerr,decoded] = rsdec(RS_codeword_test,(LT_seed_bit+LT_L+8*RS_parity_num)/8,(LT_seed_bit+LT_L)/8);
                
                decoded = double(decoded.x);
                binary_decoded = [];
                
                if (cnumerr == 0)
                    flag = 1;  % for cluster size 2, if first strand has RS ok, then skip the next strand
                else
                    binary_temp = fliplr(de2bi(decoded));
                    RS_check_indeicator = 0;
                    
                    for i = 1:(LT_seed_bit+LT_L+8*RS_parity_num)/8
                        binary_decoded = [binary_decoded binary_temp(i,:)];
                    end

                    trial_number = trial_number + 1;
                    
                    for i=1:(LT_seed_bit+LT_L+8*RS_parity_num)/8
                        for j=1:8
                            eightbit_temp(j)=LT_parity_candidate(trial_number,8*i-j+1);
                        end
                        RS_test_temp(i) = bi2de(eightbit_temp);
                    end                    

                    RS_codeword_test = gf(RS_test_temp,8);
                    [~,cnumerr,decoded] = rsdec(RS_codeword_test,(LT_seed_bit+LT_L+8*RS_parity_num)/8,(LT_seed_bit+LT_L)/8);
                    decoded = double(decoded.x);
                    
                    if (cnumerr ~= 0)
                        
                        corrected_strand = [corrected_strand ; binary_decoded];
                        binary_decoded = [];
                        
                        binary_temp = fliplr(de2bi(decoded));
                        RS_check_indeicator = 0;
                        
                        for i = 1:(LT_seed_bit+LT_L+8*RS_parity_num)/8
                            binary_decoded = [binary_decoded binary_temp(i,:)];
                        end
                        
                        corrected_strand = [corrected_strand ; binary_decoded];
                        
                        cnumerr = 5;
                        real_trial = real_trial-1;
                        
                    end
                end
                
            else
%%        
                for i=1:(LT_seed_bit+LT_L+8*RS_parity_num)/8
                    for j=1:8
                        eightbit_temp(j)=LT_parity_candidate(trial_number,8*i-j+1);
                    end
                    RS_test_temp(i) = bi2de(eightbit_temp);
                end
                RS_codeword_test = gf(RS_test_temp,8);
                [~,cnumerr,decoded] = rsdec(RS_codeword_test,(LT_seed_bit+LT_L+8*RS_parity_num)/8,(LT_seed_bit+LT_L)/8);

                decoded = double(decoded.x);
                binary_decoded = [];
            
            end
    %%        
            if(cnumerr == 0)
                RS_check_indicator = 1;
                binary_candidate = LT_parity_candidate(trial_number,:);
                real_trial = real_trial+1;

            elseif((cnumerr == 1) &&(sum(RS_test_temp ~=decoded)==1))

                binary_temp = fliplr(de2bi(decoded));
                RS_check_indeicator = 0;
                for i = 1:(LT_seed_bit+LT_L+8*RS_parity_num)/8
                    binary_decoded = [binary_decoded binary_temp(i,:)];
                end
                
                    corrected_strand = [corrected_strand ; binary_decoded];
            else
                RS_check_indeicator = 0;
                real_trial = real_trial+1;
                
            end
     %%       
            if(RS_check_indicator == 1)   
                temp_vector_candidate = fliplr(binary_candidate(1:LT_seed_bit));

                Galois_register = bi2de(temp_vector_candidate);

                rng(Galois_register);
                degree_prop = rand;
                LT_generator_temp = zeros(LT_K,1);

                for j=1:LT_K
                    if(j==1)
                        if(degree_prop < Robust_soliton(1))
                            symbolselection = datasample([1:LT_K],1,'Replace',false);
                            LT_generator_temp(symbolselection(1),1)=1;
                        end
                    end

                    if(j~=1)
                        if((sum_test(j-1)<=degree_prop) && (degree_prop < sum_test(j)))    
                            symbolselection = datasample([1:LT_K],j,'Replace',false);
                            for m=1:j
                                LT_generator_temp(symbolselection(m),1)=1;
                            end
                        end
                    end

                    if(degree_prop==1)
                        LT_generator_temp(:,1)=1;
                    end
                end

                received_xor_state = find(LT_generator_temp);  %%%find index number of 1 in LT_generator_temp

                if(sum(Inferred_idx(received_xor_state)~=LT_generator_temp(received_xor_state)) ~= 0) %% has degree other than inferred ones
                    LT_decoding_observance = LT_decoding_observance + 1;
                    subtracted_sequence = binary_candidate(LT_seed_bit+1:LT_seed_bit+LT_L);

                    for i=1:size(received_xor_state,1)
                        if(Inferred_idx(received_xor_state(i)) == 1)
                            subtracted_sequence = mod(subtracted_sequence + Inferred_bitstream(received_xor_state(i),:),2);
                            LT_generator_temp(received_xor_state(i)) = 0;
                        end              
                    end

                    test_degree_one = find(LT_generator_temp);
                    test_test_degree_size = size(test_degree_one,1);


                    if(test_test_degree_size==1)  %%% became degree 1 after inferring
                        %%%%below: we have new information from degree 1,
                        %%%%so we need a bp process

                        bp_update = 1;
                        while(bp_update == 1)
                            if(LT_real_N == 0)
                                bp_update = 0;
                                Inferred_bitstream(test_degree_one(1),:) = subtracted_sequence;
                                Inferred_idx(test_degree_one(1)) = 1;

                            end

                            if(LT_real_N ~= 0)    

                                Inferred_bitstream(test_degree_one(1),:) = subtracted_sequence;
                                Inferred_idx(test_degree_one(1)) = 1;


                                for i=1:LT_real_N
                                    if((Current_xor_number(i) ~= 0) && (Current_xor_state(i,test_degree_one(1)) == 1))
                                        Current_bitstream(i,:) = mod(Current_bitstream(i,:) + subtracted_sequence,2);
                                        Current_xor_state(i,test_degree_one(1)) = 0;
                                        Current_xor_number(i) = Current_xor_number(i) - 1;
                                    end
                                end                     
                                %%how to find the next one
                                for i=1:LT_real_N
                                    if(Current_xor_number(i) == 1)
                                        test_degree_one(1) = find(Current_xor_state(i,:));
                                        subtracted_sequence = Current_bitstream(i,:);
                                        bp_update = 1;
                                        break;                            
                                    end

                                    if(i==LT_real_N)
                                        bp_update = 0;
                                    end
                                end
                            end                               
                        end

                    else
                        LT_real_N = LT_real_N + 1;
                        Current_xor_state(LT_real_N,:) = transpose(LT_generator_temp);
                        Current_bitstream(LT_real_N,:) = subtracted_sequence;
                        Current_xor_number(LT_real_N) = test_test_degree_size;

                    end

                end

            end
            if(flag==1)
                trial_number = trial_number+1;
            end

        else
%% start strand decoding of cnumerr = 1

            corrected_size_all = size(corrected_strand);
            corrected_size = corrected_size_all(1);

            while ((k < corrected_size) && (sum(Inferred_idx)~=LT_K))
                
                k = k+1;
                real_trial =real_trial+1;
                
                binary_candidate = corrected_strand(k,:);

                temp_vector_candidate = fliplr(binary_candidate(1:LT_seed_bit));

                Galois_register = bi2de(temp_vector_candidate);

                rng(Galois_register);
                degree_prop = rand;
                LT_generator_temp = zeros(LT_K,1);


                for j=1:LT_K
                    if(j==1)
                        if(degree_prop < Robust_soliton(1))
                            symbolselection = datasample([1:LT_K],1,'Replace',false);
                            LT_generator_temp(symbolselection(1),1)=1;
                        end
                    end

                    if(j~=1)
                        if((sum_test(j-1)<=degree_prop) && (degree_prop < sum_test(j)))    
                            symbolselection = datasample([1:LT_K],j,'Replace',false);
                            for m=1:j
                                LT_generator_temp(symbolselection(m),1)=1;
                            end
                        end
                    end

                    if(degree_prop==1)
                        LT_generator_temp(:,1)=1;
                    end
                end

                received_xor_state = find(LT_generator_temp);  %%%find index number of 1 in LT_generator_temp

                if(sum(Inferred_idx(received_xor_state)~=LT_generator_temp(received_xor_state)) ~= 0)   %% has degree other than inferred ones
                    LT_decoding_observance = LT_decoding_observance + 1;
                    subtracted_sequence = binary_candidate(LT_seed_bit+1:LT_seed_bit+LT_L);

                    for i=1:size(received_xor_state,1)
                        if(Inferred_idx(received_xor_state(i)) == 1)
                            subtracted_sequence = mod(subtracted_sequence + Inferred_bitstream(received_xor_state(i),:),2);
                            LT_generator_temp(received_xor_state(i)) = 0;
                        end              
                    end

                    test_degree_one = find(LT_generator_temp);
                    test_test_degree_size = size(test_degree_one,1);


                    if(test_test_degree_size==1)  %%% became degree 1 after inferring
                        %%%%below: we have new information from degree 1,
                        %%%%so we need a bp process

                        bp_update = 1;
                        while(bp_update == 1)
                            if(LT_real_N == 0)
                                bp_update = 0;
                                Inferred_bitstream(test_degree_one(1),:) = subtracted_sequence;
                                Inferred_idx(test_degree_one(1)) = 1;

                            end

                            if(LT_real_N ~= 0)    

                                Inferred_bitstream(test_degree_one(1),:) = subtracted_sequence;
                                Inferred_idx(test_degree_one(1)) = 1;


                                for i=1:LT_real_N
                                    if((Current_xor_number(i) ~= 0) && (Current_xor_state(i,test_degree_one(1)) == 1))
                                        Current_bitstream(i,:) = mod(Current_bitstream(i,:) + subtracted_sequence,2);
                                        Current_xor_state(i,test_degree_one(1)) = 0;
                                        Current_xor_number(i) = Current_xor_number(i) - 1;
                                    end
                                end                     
                                %%how to find the next one
                                for i=1:LT_real_N
                                    if(Current_xor_number(i) == 1)
                                        test_degree_one(1) = find(Current_xor_state(i,:));
                                        subtracted_sequence = Current_bitstream(i,:);
                                        bp_update = 1;
                                        break;                            
                                    end

                                    if(i==LT_real_N)
                                        bp_update = 0;
                                    end
                                end
                            end                               
                        end

                    else
                        LT_real_N = LT_real_N + 1;
                        Current_xor_state(LT_real_N,:) = transpose(LT_generator_temp);
                        Current_bitstream(LT_real_N,:) = subtracted_sequence;
                        Current_xor_number(LT_real_N) = test_test_degree_size;

                    end

                end

            end
                
            break
                
        end
    end
    
    Trial_Number_collect(recursive-start_index+1) = trial_number;  %%save trial number
    LT_real_N_collect(recursive-start_index+1) = LT_real_N;
    LT_decoding_observance_collect(recursive-start_index+1) = LT_decoding_observance;

    if(sum(Inferred_idx)~=LT_K)
        fprintf('Final: Science_decoding_failure:inference_number_failure(0)\n')
        Result_collect(recursive-start_index+1) = 0;

    else
        if(Binary_Data_input == Inferred_bitstream)
            fprintf('Final: Science_decoding_success(1)\n')
            Result_collect(recursive-start_index+1) = 1;
        else
            fprintf('Final: Science_decoding_failure(2)\n')
            Result_collect(recursive-start_index+1) = 2;
        end   
    end
    
end

output_filename = sprintf('Decoding_final_result.txt');
[FP] = fopen(output_filename,'wt');

for i=1:end_index-start_index+1
    fprintf(FP,'%d ',Result_collect(i));
end

fclose(FP);
