with open("/tscc/nfs/home/ymadakamutil/projects/labrat/gencodecomprehensive.v28.gff3") as f:
    gff = list(f)

gff = [x for x in gff if not x.startswith('#')]

gff = [x for x in gff if 'gene_id=' in x and 'gene_name=' in x]

gff = list(map(lambda x: (x.split('gene_id=')[1].split(';')[0], x.split('gene_name=')[1].split(';')[0]), gff))

gff = dict(gff)