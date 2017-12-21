#!/usr/bin/env python3
#coding:utf-8
import os
import re
import sys

def config_fun():
    config_dict = {
        'samtools' : 'samtools',
    }
    return config_dict

# run cmd
def run_cmd(input_cmd = None):
    if input_cmd == None:
        print('Input Cmd == None!')
    print('Input Cmd:\n\t%s;' % input_cmd)
    os.system(input_cmd)

# extract head file
def extract_head(bam_file=None):
    '''
    samtools.extract_head
    :param bam_file: *.sam/.bam
    :return: head_file;
    '''
    if bam_file == None or os.path.exists(bam_file.strip()) == None:
        print('Bam File Not Find!')
        sys.exit()
    bam_file = os.path.abspath(bam_file.strip())
    head_file = re.sub('(\.bam|\.sam)$','_head.txt',bam_file)
    config_dict = config_fun()
    # cmd
    extract_cmd = ' '.join([
        config_dict['samtools'],'view','-H',bam_file,
        '>',head_file,'2>',head_file+'.log',
    ])
    run_cmd(extract_cmd)
    # return
    return head_file

# parse cigar
def cigar_parse(cigar_info=None):
    '''
    samtools.cigar_parse: parse cigar_info(45S35M4D54I243S)
    :param cigar_info: cigar_info;
    :return: ( cigar_list[ [rank_num,length,cigar_tag] ],cigar_dict{ rank_num:[length,cigar_tag], } );
    Example:
        45S34M24D
            -> list [ [45,S], [34,M], [24,D], ]
            -> dict { 0:[45,S], 1:[34,M], 2:[24,D], }
    '''
    # test
    if cigar_info == None or re.search('[MIDNSHP=X]',cigar_info) == None:
        print('Cigar Info == None');print(cigar_info)
        sys.exit()
    # parse
    cigar_info = cigar_info.strip()
    cigar_list = []
    mapping_len = '';
    for mapping_value in list(cigar_info):
        if re.search('[0-9]', mapping_value):
            mapping_len += mapping_value
        else:
            if re.search('^(M|D|I|S|H|N|X|=|P)$',mapping_value):  # match/delection/insert/soft/hard/skip_region/mismatch/match/delete_from_ref
                cigar_list.append([ mapping_len,mapping_value ])
            else:
                print('Not Find %s in CIGAR!' % mapping_value);
                print(cigar_info)
                sys.exit()
            mapping_len = ''
    cigar_dict = { rank_num : info_list for rank_num,info_list in enumerate(cigar_list) }

    return (cigar_list,cigar_dict)

if __name__ == '__main__':
    if re.match('extract_head',sys.argv[1],re.IGNORECASE):
        print( extract_head(sys.argv[2]) )
    elif re.match('cigar_parse',sys.argv[1],re.IGNORECASE):
        (cigar_list,cigar_dict) = cigar_parse(sys.argv[2])
        print(cigar_list)
        print(cigar_dict)
    else:
        print('Samtools Fun:')
        print('\tpython3 *py extract(_head) input_bam;')
        print('\tpython3 *py cigar_parse cigar_info')
    print('End')







