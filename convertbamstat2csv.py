'''
This code extracts the mapping stats from a bamstats text file <br>
it first fixes the nunber of Forward strand mapped reads <br>
Then it writes to a csv file 

'''
__author__ = "Shabhonam_caim@cegx.co.uk"
import glob


def get_records(statfile):
    '''
   
    To extract the header information the snippet code we run first.<br>
    
       
    '''

    headers = ['Total reads', 'Mapped reads', 'Forward strand', 'Reverse strand', 'Failed QC', 'Duplicates', 'Paired-end reads', "'Proper-pairs'", 'Both pairs mapped', 'Read 1', 'Read 2', 'Singletons', 'Average insert size (absolute value)', 'Median insert size (absolute value)']


    # unmapped reads(Total-mapped), Improper pairs(Both pair - proper pair), Singletons, proper-pairs:
    data ={}
    with open(statfile) as inp:
       ignore1 = inp.readline()#.strip().split(':') 
       ignore2 = inp.readline()#.strip().split(':') 
       ignore3 = inp.readline()#.strip().split(':') 
       for line in inp:
           A = line.strip().split(':') 
           if len(A) > 1:
               h, val = [A[0], A[1].split()[0]]
               assert h in headers, ' check the file should have the following headers:' + str(headers)
               data[h] = val

     
    # Fix the Forward strand values 
    data['Forward strand'] = int(data['Mapped reads']) -  int(data['Reverse strand'] )
     
    return headers, data 

    
if __name__ == '__main__':
    bwa_stat_files = glob.glob('*.bamtools_stats.txt')


    need_header = ['ID', 'Proper_pairs',  'Singletons',  'Unmapped_reads',  'Improper_pairs',  'Total_reads_ff', 'Total_reads_rr']

    pull_headers = {

        'ID' : 'ID', 
        'Proper_pairs':  'Proper_pairs', 
        'Singletons': 'Singletons' ,  
        'Unmapped_reads' : 'unmapped_reads',  
        'Improper_pairs'  : 'new_improper_pair',
        'Total_reads_ff':  'mapped_reads_ff', 
        'Total_reads_rr': 'mapped_reads_rr',
        'Forward strand:' : 'fwd_strand',
    }

    for statfile in sorted(bwa_stat_files):
        output_csv_filer = statfile[:-4] + '.csv'
        with open(output_csv_filer,'w') as outF:
            headers, data = get_records(statfile)
            outF.write(','.join(headers) + '\n')
            outF.write(','.join([str(data[h]) for h in headers]) + '\n')

    print ('DONE')     

