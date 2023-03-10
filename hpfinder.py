import sys
import gzip

def hpfinder(fq1, fq2, r1_hairpin, r2_hairpin):
    tot = 0
    found_f = 0
    found_r = 0
    with gzip.open(fq1, 'rt') as F, gzip.open(fq2, 'rt') as R:
        for (h1, s1),  (h2, s2) in zip(get_seq(F), get_seq(R)):
            tot += 1
            if r1_hairpin in s1 or r2_hairpin in s1:
               found_f +=1
            if r1_hairpin in s2 or r2_hairpin in s2:
               found_r +=1
    return tot, found_f, found_r



def get_seq(fh):
  while fh:
    header = fh.readline().strip()
    seq    = fh.readline().strip()
    sign = fh.readline()
    #qual = fh.readline()
    next(fh)
    if not header:
       break
    yield header, seq


if __name__ == '__main__':

    fq1 = sys.argv[1]
    fq2 = sys.argv[2]
    start_of_hp = sys.argv[3]

    r1_hairpin="ATGACGATGCGTTCGAGCATCGTCATT"
    r2_hairpin="AATGACGATGCTCGAACGCATCGTCAT"
    tot = 0
    found_f = 0
    found_r = 0
    with open('hairpin_dist.csv', 'w') as outF:
        outF.write('HP_length, HP_forward , HP_reverse\n')
        for num_base in range(int(start_of_hp),len(r1_hairpin)):
        
            r1_hair = r1_hairpin[:num_base] 
            r2_hair = r2_hairpin[:num_base]
            tot, found_f, found_r = hpfinder(fq1, fq2, r1_hair, r2_hair)
           
            outF.write(f'{num_base}, {found_f}, {found_r}\n')

            print (f'Length of hairpin found:{num_base}')
            print (f'found in forward :{found_f}')
            print (f'found in reverse:{found_r}')
            print ('')


        
