#!/usr/bin/env python3
__author__ = 'Haining.Wang'
#coding:utf-8
import os
import re
import sys

# Run Cmd
def run_cmd(input_cmd=None):
    '''
    :param input_cmd: Cmd Str
    :return: None
    '''
    if input_cmd == None:
        print('Input Cmd == None!')
        sys.exit()
    print('Run Cmd:\n\t%s;' % input_cmd)
    os.system(input_cmd)

# Chrom Standard Format
def chrom_format(input_chrom=None):
    '''
    :param input_chrom:chrX/Y/MT/23/24
    :return: format_chrom-X/Y/MT/1/2...
    '''
    # check
    if input_chrom == None:
        print('Input Chrom == None!')
        sys.exit()
    # strip()
    input_chrom = input_chrom.strip()
    # replace chr/CHR
    if re.match('chr',input_chrom,re.IGNORECASE):
        input_chrom = re.sub('chr','',input_chrom,re.IGNORECASE)
    # 23/34/25 => X/Y/MT
    if re.search('^23$',input_chrom):
        input_chrom = 'X'
    elif re.search('^24$',input_chrom):
        input_chrom = 'Y'
    elif re.search('^25$',input_chrom):
        input_chrom = 'MT'
    return input_chrom

# Genome Create Index
def create_index(input_genome=None):
    '''
    :param input_genome: genome;
    :return: genome.fai File;
    '''
    # check
    if input_genome != None:
        input_genome = os.path.abspath(input_genome)
    if input_genome == None or os.path.exists(input_genome) == None:
        print('Input Genome!')
        sys.exit()
    # Check Exists Fai
    if os.path.exists(input_genome+'.fai'):
        print('Have Create Index:\t%s!' % input_genome+'.fai')
        return (input_genome+'.fai')
    # Make Index
    else:
        index_cmd = ' '.join([
            'samtools', 'faidx', input_genome,
            '>', input_genome+'.fai.log',
        ])
        run_cmd(index_cmd)
        return input_genome+'.fai'

# Read Fai to fai_dict
def read_fai(input_fai=None):
    '''
    :param input_fai:fai_file;
    :return: fai_dict:{ chrom:[ chr_total_len, chrom_star, line_len, line_len_plus_space ], ...,};
    '''
    # ref   ref_len start   line_length line_length_plus_\n
    # NW_003613580.1  8779783 116     80      81
    # NW_003613581.1  8081566 8889763 80      81
    # NW_003613582.1  6666273 17072465        80      81
    # check
    if input_fai == None:
        print('Input Fai == None!')
        sys.exit()
    # Get Fai Dict
    fai_dict = {}
    with open(input_fai,'r') as fai_h:
        for line in fai_h:
            infos = re.split('\t',line.strip())
            infos[0] = chrom_format(infos[0]) # format chrom
            if len(infos) != 5: # check fai is Correct
                print('Check Fai File, line_len must > 0!')
                sys.exit()
            fai_dict[ infos[0] ] = list(map( lambda x:int(x),infos[1:] ))
    fai_h.close()

    return fai_dict

# Extract base from pos and genome
def extract_base(input_pos=None,input_genome=None,plus_len=100):
    '''
    :param input_pos: chr1:pos1(:pos2), Note: sep-':';
    :param input_genome: genome_file;
    :param plus_len: int-plus_len...target_bases....plus_len;
    :return: base_file;
    '''
    # check
    if input_pos == None or input_genome == None:
        print('Input Pos or Input_Genome == None!')
        sys.exit()
    # check input_pos
    pos_list = re.split( ':',re.sub(',','',input_pos) )
    if len(pos_list) < 2:
        print('Input_Pos Format: chrom:pos1(:pos2)!')
        sys.exit()
    for pos_num,pos_value in enumerate(pos_list):
        if pos_num >= 1 and int(pos_value) < 0:
                print('Input Pos Only >= 0!')
                sys.exit()
    # check and read fai
    input_genome = os.path.abspath(input_genome)
    fai_file = input_genome+'.fai'
    if os.path.exists(fai_file) == None:
        create_index(input_genome)
    fai_dict = read_fai(fai_file)

    # extract fa
    genome_h = open(input_genome,'r')
    pos_list[0] = chrom_format(pos_list[0])
    # chrom not in fai
    if pos_list[0] not in fai_dict:
        print('Chrom-%s Not In Fai File-%s!' % (pos_list,fai_file) )
        sys.exit()
    else:
        # fai info
        (chrom_len,chrom_start,line_len,line_space_len) = list(map( lambda x:int(x),fai_dict[pos_list[0]] ))
        line_space_len = line_space_len-line_len
        chrom_end = chrom_start + chrom_len + chrom_len//line_len
        extract_pos = sorted( list(map(lambda x: int(x), pos_list[1:])) ); #print(extract_pos);sys.exit()
        # only one pos => start...end
        if len(extract_pos) == 1:
            extract_pos.append(extract_pos[0])
        # pos <= chrom_len and remap pos to genome file
        for pos_num,pos_value in enumerate(extract_pos):
            if pos_value > chrom_len: # check pos < chrom_total_len
                print('Input Pos > Chrom_total_len');
                pos_value = chrom_len
            space_len = pos_value//line_len if pos_value%line_len != 0 else (pos_value//line_len) -1
            extract_pos[pos_num] = pos_value+chrom_start+space_len*line_space_len
        move_pos = extract_pos[0]-plus_len if extract_pos[0]-plus_len > chrom_start else chrom_start # Note chrom_start
        # output
        output_fold = os.getcwd()
        output_file_name = os.path.basename(input_genome)+'.%s.fa' \
            % '_'.join(pos_list)
        output_file = os.path.join( output_fold,output_file_name )
        output_h = open(output_file,'w')
        genome_h.seek(move_pos,0)
        for line in genome_h:
            for base_value in list(line):
                if move_pos+1 == extract_pos[0]:
                    output_h.write('*')
                if len(base_value.strip()) == 1:
                    output_h.write(base_value)
                if move_pos+1 == extract_pos[1]:
                    output_h.write('*')
                move_pos += 1
                if move_pos > extract_pos[1]+plus_len or move_pos > chrom_end:
                    break
            if move_pos > extract_pos[1]+plus_len or move_pos > chrom_end:
                break
        output_h.close()
        genome_h.close()

        return output_file

if __name__ == '__main__':
    if re.search('create_index',sys.argv[1],re.IGNORECASE):
        print( create_index(sys.argv[2]) )
    elif re.search('extract_base',sys.argv[1],re.IGNORECASE):
        if len(sys.argv) == 4:
            sys.argv.append(100)
        print( extract_base(sys.argv[2],sys.argv[3],sys.argv[4]) )
    else:
        print('Fun:')
        print('\tcreate_index genome_file;')
        print('\textract_base chr1:pos1(:pos2) (100);')
        sys.exit()
    print('End')


