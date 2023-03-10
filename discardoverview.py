import glob


def enigma_loader(statfile, direction ='ff'):
    if direction == 'ff':
       enigmafile = statfile.replace('bwameth.bamtools_stats.txt', 'bwameth.enigmaff.bamtools_stats.txt') 
    else:
       enigmafile = statfile.replace('bwameth.bamtools_stats.txt', 'bwameth.enigmarr.bamtools_stats.txt') 

    data ={}
    with open(enigmafile) as inp:
       ignore1 = inp.readline()#.strip().split(':') 
       ignore2 = inp.readline()#.strip().split(':') 
       ignore3 = inp.readline()#.strip().split(':') 
       for line in inp:
            
           A = line.strip().split(':')
           if len(A) > 1:
             
             head, tail = A[0], A[1]
             if head in 'Mapped reads:':
                data['mapped_reads'] = int(tail.split()[0]) 

    return data


def get_records(statfile):
    '''
    Total reads:       5920478
    Mapped reads:      5821730      (98.3321%)
    Forward strand:    2998601      (50.648%)
    Reverse strand:    2921877      (49.352%)
    Failed QC:         220790       (3.72926%)
    Duplicates:        0    (0%)
    Paired-end reads:  5920478      (100%)
    'Proper-pairs':    5217601      (88.128%)
    Both pairs mapped: 5799462      (97.956%)
    Read 1:            2959462
    Read 2:            2961016
    Singletons:        22268        (0.376118%)
    Average insert size (absolute value): 21488.2
    Median insert size (absolute value): 210
    unmapped reads(Total-mapped), Improper pairs(Both pair - proper pair), Singletons, proper-pairs:

    '''
    data ={}
    with open(statfile) as inp:
       ignore1 = inp.readline()#.strip().split(':') 
       ignore2 = inp.readline()#.strip().split(':') 
       ignore3 = inp.readline()#.strip().split(':') 
       for line in inp:
            
         A = line.strip().split(':')
         if len(A) >  1:
           head, tail = A[0], A[1]
           if head in 'Total reads:':
              data['total_reads'] = int(tail.split()[0]) 
           if head in "'Proper-pairs'":
              data['Proper_pairs'] = int(tail.split()[0]) 
           if head in 'Singletons:':
              data['Singletons']  = int(tail.split()[0]) 

           if head in 'Mapped reads:':
              data['mapped_reads'] = int(tail.split()[0]) 

           if head in 'Reverse strand:':
              data['rev_strand'] = int(tail.split()[0]) 
   
           if head in 'Both pairs mapped:':
              data['both_mapped'] = int(tail.split()[0]) 
    #-----
    data['correct_forward'] = data['mapped_reads'] -  data['rev_strand'] 
    data['unmapped_reads']  = data['total_reads']  -  data['mapped_reads'] 
    data['improper_pair']   = data['both_mapped'] -  data['Proper_pairs']
    return data 

    

bwa_stat_files = glob.glob('*_discarded_R1.bwameth.bamtools_stats.txt')


need_header = ['ID', 'Proper_pairs',  'Singletons',  'Unmapped_reads',  'Improper_pairs',  'Total_reads_ff', 'Total_reads_rr']

pull_headers = {

    'ID' : 'ID', 
    'Proper_pairs':  'Proper_pairs', 
    'Singletons': 'Singletons' ,  
    'Unmapped_reads' : 'unmapped_reads',  
    'Improper_pairs'  : 'new_improper_pair',
    'Total_reads_ff':  'mapped_reads_ff', 
     'Total_reads_rr': 'mapped_reads_rr',
}


with open('discarded_file_overview.csv','w') as outF:
 pipe =  (','.join( [ h for h in need_header] ))
 outF.write(pipe + '\n')
    
 for statfile in sorted(bwa_stat_files):
   
    data = get_records(statfile)
    data['ID'] = statfile.replace('_discarded_R1.bwameth.bamtools_stats.txt','')

    try:
      data['mapped_reads_ff'] = enigma_loader(statfile, direction ='ff')['mapped_reads']
      data['mapped_reads_rr'] = enigma_loader(statfile, direction ='rr')['mapped_reads']
    except:
      print('Enigma_file_missing for sample:', data['ID'])
      continue
      
    data['new_improper_pair']   =  data['improper_pair'] -  (data['mapped_reads_ff'] +  data['mapped_reads_rr'])

    sample_name  = statfile.replace('bwameth.bamtools_stats.txt', 'bwameth.enigmaff.bamtools_stats.txt') 


    pipe =  (','.join( [ str(data[pull_headers[h]]) for h in need_header] ))
    outF.write(pipe + '\n')


    