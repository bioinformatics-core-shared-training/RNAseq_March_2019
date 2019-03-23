import sys
import csv

def main(gffcompare_annotated_gtf, gffcompare_tracking ):

  for row in csv.reader(open(gffcompare_annotated_gtf), delimiter="\t"):
    print(row)
  
  


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
