#!/usr/bin/env python
#!/bin/sh


import pandas as pd

def convert_to_gff(bed_file, out_file):
    Bed= pd.read_csv(bed_file, sep= '\t', header= -1) 
    name= Bed[0]
    start= Bed[1]
    end= Bed[2]
    geneId= ' '
    score= '1'
    out = open(out_file, "w")
    for i in range(len(Bed.index)):
        StartinGFF= start[i]
        EndinGFF= end[i]
        NameinGFF= name[i]
        ScoreinGFF= score
        AttrinGFF= geneId
        if len(Bed.columns) > 3:
            geneId= Bed[3]
            AttrinGFF= geneId[i]
            score= Bed[4] 
            ScoreinGFF= score[i]
        Format= ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (NameinGFF, 'MACE', ' ', StartinGFF, EndinGFF, ScoreinGFF, ' ', ' ', AttrinGFF))
        out.write(Format+'\n')
    out.close()
    
if __name__ == "__main__":
    import pandas as pd
    from argparse import ArgumentParser
    from os.path import split
    parser = ArgumentParser("convert a Bed to a gff file")
    parser.add_argument("bed_file")
    parser.add_argument("out_file")
    args = parser.parse_args()
    convert_to_gff(args.bed_file, args.out_file)
