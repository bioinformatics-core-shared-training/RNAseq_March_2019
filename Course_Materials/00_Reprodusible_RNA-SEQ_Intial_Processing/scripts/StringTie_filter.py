import sys
import csv

def main(gffcompare_annotated_gtf, gffcompare_tracking ):

  for row in csv.reader(open(gffcompare_tracking), delimiter="\t"):
    
    stringtie_class = row[3]
    transcript = row[4].split("|")[1]

    print(transcript, stringtie_class)
  
#  for row in csv.reader(open(gffcompare_annotated_gtf), delimiter="\t"):
#    print(row)
  
  


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
