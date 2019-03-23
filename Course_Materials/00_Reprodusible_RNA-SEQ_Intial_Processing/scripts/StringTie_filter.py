import sys
import csv
from collections import defaultdict

def main(gffcompare_annotated_gtf, gffcompare_tracking, exon_number_filter, class_filter ):

  transcript_class_dict = {}
  tag_dict = defaultdict(dict)
  
  for row in csv.reader(open(gffcompare_tracking), delimiter="\t"):
    
    stringtie_class = row[3]
    transcript = row[4].split("|")[1]

    transcript_class_dict[transcript] = stringtie_class
    
  
  for row in csv.reader(open(gffcompare_annotated_gtf), delimiter="\t"):


    for tag in row[-1].split(";"):
      tag_value = tag.strip(" ").split(" ")
      
      transcript_id = ""
      
      if len(tag_value)==2:
        tag, value = tag_value
        
        print(tag, value)
        
        if tag=="transcript_id":
          transcript_id=value
        else:
          tag_dict[transcript_id][tag] = value
        
    #print(tag_dict)

      
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
