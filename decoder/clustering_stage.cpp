/*
* Copyright (c) 2021 Jaeho Jeong from Coding and Cryptography Lab (CCL),
* Department of Electrical and Computer Engineering, Seoul National University, South Korea.
* Email address: jaehoj@ccl.snu.ac.kr
* All Rights Reserved.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <math.h>
#define LT_N 18000
#define LEN 152
#define MIN_DISTANCE 2
typedef struct oligo {
	char *s;
	int oc_index;	//oligo_cluster_index_number
	int cnt;
	double avg_mul_qscore;
}oligo;

bool sort_by_oc(oligo a, oligo b)
{
	return (a.oc_index < b.oc_index);
}
bool sort_by_cnt(oligo a, oligo b)
{
	return (a.cnt > b.cnt);
}
bool sort_by_qscore_mul(oligo a, oligo b)
{
	return (a.avg_mul_qscore > b.avg_mul_qscore);
}

int hamming_dist(char *input1, char *input2, int length)
{
	int cnt = 0;
	for (int i = 0; i < length; i++)
	{
		if (input1[i] != input2[i]) cnt++;
	}
	return cnt;
}

int max(int a, int b, int c, int d)
{
	int x, y;
	if (a > b) x = a;
	else x = b;
	if (c > d) y = c;
	else y = d;
	if (x > y) return x;
	else return y;
}

int main(int argc, char *argv[])
{
	int pos = 1;
	char *line1_1, *line1_2, *line2, *line3, *line4;
	double *Qscore;
	double mul_qscore;
	int i, j, k, loop_cnt, line_num, N_count, ham, oc_index_num, index_top, total_oc, target_oc_index, total_discard;
	int num_A, num_C, num_G, num_T, flag, temp_num, max_num, location_one, k1, k2, k3, k4, flag_for_cluster2;

	FILE *in1, *in2;
	FILE *out, *out1, *out2, *out_discard, *out_cnt;

	char buf[100];
	int index_number;
	
	index_number = atoi(argv[pos++]);
	k1 = atoi(argv[pos++]);
	k2 = atoi(argv[pos++]);
	k3 = atoi(argv[pos++]);
	k4 = atoi(argv[pos++]);
	printf("Input file name? ");
	scanf("%s", buf);
	in1 = fopen(buf, "r");
	in2 = fopen(buf, "r");


	sprintf(buf, "raw_read%d.txt", index_number);
	out = fopen(buf, "w");
	sprintf(buf, "raw_read_symbol_cnt%d.txt", index_number);
	out1 = fopen(buf, "w");
	sprintf(buf, "raw_read_symbol_cnt_ocNum%d.txt", index_number);
	out2 = fopen(buf, "w");
	sprintf(buf, "raw_read_symbol_cnt_discarded%d.txt", index_number);
	out_discard = fopen(buf, "w");
	sprintf(buf, "raw_cnt%d.txt", index_number);
	out_cnt = fopen(buf, "w");

	////////////////////////////
	line1_1 = (char*)malloc(sizeof(char) * (LEN + 2));
	line1_2 = (char*)malloc(sizeof(char) * (LEN + 2));
	line2 = (char*)malloc(sizeof(char) * (LEN + 2));
	line3 = (char*)malloc(sizeof(char) * (LEN + 2));
	line4 = (char*)malloc(sizeof(char) * (LEN + 2));

	//////////////////////////// calculate q score
	Qscore = (double*)malloc(sizeof(double) * 41);
	for (i = 0; i < 41; i++)
	{
		Qscore[i] = 1.0 - pow(10.0, (double)(i) / (-10.0));
	}

	////////////////////////////check how many lines in FASTQ file
	line_num = 0;
	while (fscanf(in1, "%s", line1_1) != EOF)
	{
		fscanf(in1, "%s", line1_2);
		fscanf(in1, "%s", line2);
		fscanf(in1, "%s", line3);
		fscanf(in1, "%s", line4);
		if (line_num % 10000 == 0) printf("checking how many lines in fastq file: %d\n", line_num);
		line_num++;
	}
	std::fclose(in1);
	////////////////////////////
	oligo *read;
	read = (oligo*)malloc(sizeof(oligo)*line_num);
	int *cnt_array, *index_array;
	cnt_array = (int*)malloc(sizeof(int)*line_num);
	index_array = (int*)malloc(sizeof(int)*line_num);

	for (i = 0; i < line_num; i++)
	{
		read[i].s = (char*)malloc(sizeof(char)*(LEN + 2));
		cnt_array[i] = 0;
	}


	////////////////////////////

	loop_cnt = 0;
	N_count = 0;
	oc_index_num = 0;

	printf("line number: %d\n\n", line_num);
	while (fscanf(in2, "%s", line1_1) != EOF)
	{
		fscanf(in2, "%s", line1_2);
		fscanf(in2, "%s", line2);
		fscanf(in2, "%s", line3);
		fscanf(in2, "%s", line4);

		////////////////if N found, then skip
		for (i = 0; i < LEN; i++)
		{
			if (line2[i] == 'N')
			{
				N_count++;
				break;
			}
		}
		if (i != LEN)
		{
			line_num--;
			continue;
		}

		////////////////line1, line2, line3
		mul_qscore = 1.0;
		for (i = 0; i < LEN; i++)
		{
			mul_qscore *= Qscore[(int)line4[i] - 33];
		}

		strncpy(read[loop_cnt].s, line2, LEN);
		read[loop_cnt].avg_mul_qscore = mul_qscore;

		index_top = 0;
		for (i = 0; i < loop_cnt; i++)
		{
			ham = hamming_dist(read[loop_cnt].s, read[i].s, LEN);
			if (ham <= MIN_DISTANCE)
			{
				index_array[index_top] = i;
				index_top++;
			}
		}

		if (index_top > 0)	//has similar distacnce strands -> cluster them together (C1)
		{
			read[loop_cnt].oc_index = read[index_array[0]].oc_index;
			cnt_array[read[loop_cnt].oc_index]++;
			for (i = 1; i < index_top; i++)
			{
				if (read[index_array[i]].oc_index != read[loop_cnt].oc_index)
				{
					target_oc_index = read[index_array[i]].oc_index;
					cnt_array[read[loop_cnt].oc_index] += cnt_array[target_oc_index];
					cnt_array[target_oc_index] = 0;
					for (j = 0; j < loop_cnt; j++)
					{
						if (read[j].oc_index == target_oc_index)
						{
							read[j].oc_index = read[loop_cnt].oc_index;
						}
					}
				}
			}
		}
		else//new strands
		{
			read[loop_cnt].oc_index = oc_index_num;
			cnt_array[oc_index_num]++;
			oc_index_num++;
		}
		//////////////////line3, line4

		loop_cnt++;
		if (loop_cnt % 1000 == 0) printf("loop number: %d / %d\n", loop_cnt, line_num);
	}
	std::fclose(in2);

	std::stable_sort(read, read + loop_cnt, sort_by_oc);

	for (i = 0; i < loop_cnt; i++)
	{
		read[i].cnt = cnt_array[read[i].oc_index];
	}
	std::stable_sort(read, read + loop_cnt, sort_by_cnt);

	printf("line number without N: %d\n\n", line_num);

	for (i = 0; i < loop_cnt; i++)
	{
		if (read[i].cnt == 1)
		{
			location_one = i;
			break;
		}
	}
	std::stable_sort(read + location_one, read + loop_cnt, sort_by_qscore_mul); //sort by qscore (C3)

	////////////////////////////////// cluster size >= 3 then majority decision
	char* print_strand;	//save consensus at print_strand
	print_strand = (char*)malloc(sizeof(char)*LEN);
	total_oc = 0;
	total_discard = 0;
	
	int min_oligo_dist, oligo_dist;

	for (i = 0; i < line_num; i)
	{
		if (i % 10000 == 0) printf("clustering line number: %d / %d\n", i, line_num);
		flag = 0;
		for (j = 0; j < LEN; j++)
		{
			num_A = 0;
			num_C = 0;
			num_G = 0;
			num_T = 0;
			temp_num = 0;
			for (k = 0; k < read[i].cnt; k++)
			{
				if (read[i + k].s[j] == 'A') num_A++;
				else if (read[i + k].s[j] == 'C') num_C++;
				else if (read[i + k].s[j] == 'G') num_G++;
				else if (read[i + k].s[j] == 'T') num_T++;
			}
			max_num = max(num_A, num_C, num_G, num_T);
			if (max_num == num_A)
			{
				temp_num++;
				print_strand[j] = 'A';
			}
			if (max_num == num_C)
			{
				temp_num++;
				print_strand[j] = 'C';
			}
			if (max_num == num_G)
			{
				temp_num++;
				print_strand[j] = 'G';
			}
			if (max_num == num_T)
			{
				temp_num++;
				print_strand[j] = 'T';
			}
			if (temp_num != 1) flag++;
			else
			{
				for (k = 0; k < read[i].cnt; k++)
				{
					read[i + k].s[j] = print_strand[j];
				}
			}
		}
		if (flag == 0) //max likelihood
		{
			//////////////discard reads with k1<=oligo_dist<=k2 or k3<=oligo_dist<=k4 (C2)
			if (read[i].cnt == 1)
			{
				min_oligo_dist = LEN;
				for (int ii = 0; ii < location_one; ii)
				{
					oligo_dist = hamming_dist(read[i].s, read[ii].s, LEN);
					if (min_oligo_dist > oligo_dist) min_oligo_dist = oligo_dist;
					//ii += read[ii].cnt;
					ii++;
				}
				//if (k1 <= min_oligo_dist && min_oligo_dist <= k2)
				if ((k1 <= min_oligo_dist && min_oligo_dist <= k2) || (k3 <= min_oligo_dist && min_oligo_dist <= k4))
				{
					fprintf(out_discard, "%152.152s\t%d\t%d\n", read[i].s, read[i].cnt, min_oligo_dist);
					i++;
					total_discard++;
					continue;
				}
			}
			/////////////////////////////////////////
			for (j = 0; j < LEN; j++)
			{
				//save consensus at print_strand
				switch (print_strand[j])
				{
				case 'A':
					fprintf(out, "0 0 ");
					fprintf(out1, "A");
					break;
				case 'C':
					fprintf(out, "0 1 ");
					fprintf(out1, "C");
					break;
				case 'G':
					fprintf(out, "1 0 ");
					fprintf(out1, "G");
					break;
				case 'T':
					fprintf(out, "1 1 ");
					fprintf(out1, "T");
					break;
				case 'N':
					fprintf(out, "N N ");//should not exist
					fprintf(out1, "N");//should not exist
					break;
				default:
					break;
				}
			}
			fprintf(out, "\n");
			fprintf(out1, "\t%d\n", read[i].cnt);
			fprintf(out_cnt, "0\n");
		}
		else
		{
			flag_for_cluster2 = 1;
			for (k = 0; k < read[i].cnt; k++)
			{
				for (j = 0; j < LEN; j++)
				{
					switch (read[i + k].s[j])
					{
					case 'A':
						fprintf(out, "0 0 ");
						fprintf(out1, "A");
						break;
					case 'C':
						fprintf(out, "0 1 ");
						fprintf(out1, "C");
						break;
					case 'G':
						fprintf(out, "1 0 ");
						fprintf(out1, "G");
						break;
					case 'T':
						fprintf(out, "1 1 ");
						fprintf(out1, "T");
						break;
					case 'N':
						fprintf(out, "N N ");//should not exist
						fprintf(out1, "N");//should not exist
						break;
					default:
						break;
					}
				}
				fprintf(out, "\n");
				fprintf(out1, "\t%d\toverlapped\n", read[i + k].cnt);
				if (read[i + k].cnt == 2)
				{
					fprintf(out_cnt, "%d\n", flag_for_cluster2);
					flag_for_cluster2++;
				}
				else fprintf(out_cnt, "0\n");
			}
		}
		i += read[i].cnt;
		total_oc++;
	}

	//////////////////////////////////
	for (i = 0; i < loop_cnt; i++) fprintf(out2, "%152.152s\t%d\t%d\t%.5f\n", read[i].s, read[i].oc_index, read[i].cnt, read[i].avg_mul_qscore);

	fprintf(out1, "line_number: %d\nline number without N: %d\nstrand number: %d\ndiscard number: %d\n", line_num + N_count, line_num, total_oc, total_discard);
	fprintf(out_discard, "line_number: %d\nline number without N: %d\nstrand number: %d\ndiscard number: %d\n", line_num + N_count, line_num, total_oc, total_discard);

	printf("\nline number: %d\n", line_num + N_count);
	printf("line number without N: %d\n", line_num);
	printf("strand number: %d\n", total_oc);
	printf("discard number: %d\n", total_discard);

	std::fclose(out);
	std::fclose(out1);
	std::fclose(out2);
	std::fclose(out_discard);
	std::fclose(out_cnt);
	return 0;
}