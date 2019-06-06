from mygene import MyGeneInfo
from pprint import pprint
import csv
import sys

mg = MyGeneInfo()
dict_symbol = {}

def gene_name(ensg_id):
	gene = mg.getgene(ensg_id,fields='symbol')
	if gene != None and type(gene) is dict:
		return gene['symbol']
	elif type(gene) is list:
		print(ensg_id)
		pprint(gene)
		return ensg_id
	else:
		return ensg_id

'''
print(gene_name('ENSG00000273842'))
'''

with open('sample.txt','r') as f:
	rows = csv.reader(f, delimiter='\t')
	with open('gene_table.csv','w') as f_write:
		f_csv = csv.writer(f_write)
		i = 0
		for row in rows:
			i += 1
			name = gene_name(row[0][:15])
			f_csv.writerow([row[0][:15],name])
			sys.stdout.write(str(i)+'/60483\r')
			sys.stdout.flush()

