import sys
import csv

def main(gffcompare_annotated_gtf, gffcompare_tracking ):

  transcript_class_dict = {}
  
  for row in csv.reader(open(gffcompare_tracking), delimiter="\t"):
    
    stringtie_class = row[3]
    transcript = row[4].split("|")[1]

    transcript_class_dict[transcript] = stringtie_class

  
  for row in csv.reader(open(gffcompare_annotated_gtf), delimiter="\t"):
    
    print(row)
    
    for tag in row[-1].split(";"):
      print(tag)
  
  


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
