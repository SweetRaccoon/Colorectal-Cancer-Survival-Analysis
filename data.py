import json
import os
import csv
from pprint import pprint
from xml.etree import ElementTree


# class definition 

class Patient_data(object):
	pat_count = 0;

	def __init__(self,patient_id):
		self.patient_id = patient_id
		self.tp53_expression = None
		self.days_to_last_followup = []
		self.days_to_death = []
		self.days_to_new_tumor_event_after_initial_treatment = []
		self.vital_status = []
		self.days_to_last_known_alive = []
		self.tp53_mutation = []
		self.valid_varient_classification = None
		Patient_data.pat_count += 1

	def display_count(self):
		print(Patient_data.pat_count)

	def display_data(self):
		pprint(self.__dict__)

	def retrieve_tp53_expression(self,my_dict):
		if self.patient_id in my_dict:
			with open('RNA-Seq/Exp_Data/'+my_dict[self.patient_id],'r') as f:
				pat_dict = csv.reader(f, delimiter='\t')
				for row in pat_dict:
					if row[0].startswith('ENSG00000141510'):
						self.tp53_expression = row[1]
		if self.tp53_expression == None:
			self.tp53_expression = 'unknown'

	def fix_mutation(self):
		if len(self.tp53_mutation) == 0:
			self.tp53_mutation.append('unknown')



# data definition

dict_patient = {}


namespaces = {'coad':'http://tcga.nci/bcr/xml/clinical/coad/2.7',
	'xsi':'http://www.w3.org/2001/XMLSchema-instance',
	'admin':'http://tcga.nci/bcr/xml/administration/2.7',
	'clin_shared':'http://tcga.nci/bcr/xml/clinical/shared/2.7',
	'shared':'http://tcga.nci/bcr/xml/shared/2.7',
	'shared_stage':'http://tcga.nci/bcr/xml/clinical/shared/stage/2.7',
	'coad_read_shared':'http://tcga.nci/bcr/xml/clinical/shared/coad_read/2.7' ,
	'coad_nte':'http://tcga.nci/bcr/xml/clinical/coad/shared/new_tumor_event/2.7/1.0',
	'nte':'http://tcga.nci/bcr/xml/clinical/shared/new_tumor_event/2.7',
	'rx':'http://tcga.nci/bcr/xml/clinical/pharmaceutical/2.7',
	'rad':'http://tcga.nci/bcr/xml/clinical/radiation/2.7',
	'follow_up_v1.0':'http://tcga.nci/bcr/xml/clinical/coad/followup/2.7/1.0'}

# clinical

path = 'clinical/data'
files = os.listdir(path)
for fn in files:
	pat_id = fn[41:45]
	dict_patient[pat_id] = Patient_data(pat_id)
	with open(path+'/'+fn,'r') as f:
		tree = ElementTree.parse(f)
		root = tree.getroot()
		patient = root.find('coad:patient', namespaces)

		# deal with initials

		vital_status = patient.find('clin_shared:vital_status', namespaces)
		dict_patient[pat_id].vital_status.append(vital_status.text)

		days_to_death = patient.find('clin_shared:days_to_death', namespaces)
		if days_to_death.text == None:
			dict_patient[pat_id].days_to_death.append('unknown')
		else:
			dict_patient[pat_id].days_to_death.append(days_to_death.text)

		days_to_last_known_alive = patient.find('clin_shared:days_to_last_known_alive', namespaces)
		if days_to_last_known_alive.text == None:
			dict_patient[pat_id].days_to_last_known_alive.append('unknown')
		else:
			dict_patient[pat_id].days_to_last_known_alive.append(days_to_last_known_alive.text)

		days_to_last_followup = patient.find('clin_shared:days_to_last_followup', namespaces)
		if days_to_last_followup.text == None:
			dict_patient[pat_id].days_to_last_followup.append('unknown')
		else:
			dict_patient[pat_id].days_to_last_followup.append(days_to_last_followup.text)

		new_tumor_events = patient.find('coad_nte:new_tumor_events', namespaces)
		tumor_event_flag = new_tumor_events.find('nte:new_tumor_event_after_initial_treatment', namespaces)
		if tumor_event_flag.text == 'YES':
			for each in new_tumor_events.findall('coad_nte:new_tumor_event', namespaces):
				days_to_new_tumor_event_after_initial_treatment = each.find('nte:days_to_new_tumor_event_after_initial_treatment', namespaces)
				dict_patient[pat_id].days_to_new_tumor_event_after_initial_treatment.append(days_to_new_tumor_event_after_initial_treatment.text)
		else:
			dict_patient[pat_id].days_to_new_tumor_event_after_initial_treatment.append('unknown')

		# deal with follow ups

		flups = patient.find('coad:follow_ups', namespaces)
		for flup in flups.iterfind('follow_up_v1.0:follow_up', namespaces):
			vital_status = flup.find('clin_shared:vital_status', namespaces)
			dict_patient[pat_id].vital_status.append(vital_status.text)

			days_to_death = flup.find('clin_shared:days_to_death', namespaces)
			if days_to_death.text == None:
				dict_patient[pat_id].days_to_death.append('unknown')
			else:
				dict_patient[pat_id].days_to_death.append(days_to_death.text)

			days_to_last_known_alive = flup.find('clin_shared:days_to_last_known_alive', namespaces)
			if days_to_last_known_alive != None:
				if days_to_last_known_alive.text == None:
					dict_patient[pat_id].days_to_last_known_alive.append('unknown')
				else:
					dict_patient[pat_id].days_to_last_known_alive.append(days_to_last_known_alive.text)

			days_to_last_followup = flup.find('clin_shared:days_to_last_followup', namespaces)
			if days_to_last_followup.text == None:
				dict_patient[pat_id].days_to_last_followup.append('unknown')
			else:
				dict_patient[pat_id].days_to_last_followup.append(days_to_last_followup.text)

			new_tumor_events = flup.find('coad_nte:new_tumor_events', namespaces)
			tumor_event_flag = new_tumor_events.find('nte:new_tumor_event_after_initial_treatment', namespaces)
			if tumor_event_flag.text == 'YES':
				for each in new_tumor_events.findall('coad_nte:new_tumor_event', namespaces):
					days_to_new_tumor_event_after_initial_treatment = each.find('nte:days_to_new_tumor_event_after_initial_treatment', namespaces)
					dict_patient[pat_id].days_to_new_tumor_event_after_initial_treatment.append(days_to_new_tumor_event_after_initial_treatment.text)

# tp53 expression

dict_filename = {}

# open metadata and write filename to dict

with open('RNA-Seq/metadata.json','r') as f_data:
	data = json.load(f_data)

	for index in range(len(data)):
		cdata = data[index]
		fn = cdata['file_name']
		if fn.endswith('FPKM-UQ.txt.gz'):
			fn = fn[0:48]
			pat_id = cdata['cases'][0]['submitter_id'][8:12]
			if pat_id in dict_patient:
				dict_filename[pat_id] = fn

# retrieve expression

for key in dict_patient:
	dict_patient[key].retrieve_tp53_expression(dict_filename)
	print('ID:',dict_patient[key].patient_id,'tp53 expression retrieved successful')

# clear dict

dict_filename.clear()




# retrieve mutation

filepath = 'Simple_Nucleotide_Variation/gdc_download_20170531_155504/TCGA.COAD.mutect.af65d530-7976-4cd0-8ec5-2af0f4dbb3a6.DR-6.0.somatic.maf'

with open(filepath, 'r') as f:
	mutation_table = csv.reader(f, delimiter='\t')
	for x in range(0,4):
		next(mutation_table)
	headers = next(mutation_table)
	for row in mutation_table:
		if row[0] == 'TP53':
			pat_id = row[15][8:12]
			if pat_id in dict_patient:
				dict_patient[pat_id].tp53_mutation.append(row[8])

for key in dict_patient:
	dict_patient[key].fix_mutation()

# get max and min and valid var and max d_death and last vital status

temp_list = []

valid_vars = 	['Frame_Shift_Del', 
				'Frame_Shift_Ins',
				'In_Frame_Del', 
				'In_Frame_Ins', 
				'Missense_Mutation', 
				'Nonsense_Mutation', 
				'Silent', 
				'Splice_Site', 
				'Translation_Start_Site', 
				'Nonstop_Mutation']

for key in dict_patient:
	temp_list = []
	for index in range(len(dict_patient[key].days_to_last_followup)):
		if dict_patient[key].days_to_last_followup[index] != 'unknown':
			temp_list.append(int(dict_patient[key].days_to_last_followup[index]))
	if len(temp_list) != 0:
		dict_patient[key].max_followup = max(temp_list)
	else:
		dict_patient[key].max_followup = 'unknown'

	temp_list = []
	for index in range(len(dict_patient[key].days_to_new_tumor_event_after_initial_treatment)):
		if dict_patient[key].days_to_new_tumor_event_after_initial_treatment[index] != 'unknown':
			temp_list.append(int(dict_patient[key].days_to_new_tumor_event_after_initial_treatment[index]))
	if len(temp_list) != 0:
		dict_patient[key].min_new_tumor = min(temp_list)
	else:
		dict_patient[key].min_new_tumor = 'unknown'
	for index in range(len(dict_patient[key].tp53_mutation)):
		if dict_patient[key].tp53_mutation[index] in valid_vars:
			dict_patient[key].valid_varient_classification = 'Yes'
	if dict_patient[key].valid_varient_classification == None:
		dict_patient[key].valid_varient_classification = 'No'

	temp_list = []
	for index in range(len(dict_patient[key].days_to_death)):
		if dict_patient[key].days_to_death[index] != 'unknown':
			temp_list.append(int(dict_patient[key].days_to_death[index]))
	if len(temp_list) != 0:
		dict_patient[key].max_death = max(temp_list)
	else:
		dict_patient[key].max_death = 'unknown'
	# following code may have bug with all Nones
	i = -1
	temp = dict_patient[key].vital_status[i]
	while temp == None :
		i -= 1
		temp = dict_patient[key].vital_status[i]
	dict_patient[key].last_vital_status = temp


	if temp == 'Dead':
		if dict_patient[key].max_death != 'unknown':
			dict_patient[key].overall_survival = dict_patient[key].max_death 
		else:
			dict_patient[key].overall_survival = 'unknown'
	else:
		if dict_patient[key].days_to_last_known_alive[-1] != 'unknown':
			dict_patient[key].overall_survival = dict_patient[key].days_to_last_known_alive[-1]
		else:
			dict_patient[key].overall_survival = dict_patient[key].max_followup

	if temp == 'Dead':
		if dict_patient[key].min_new_tumor != 'unknown': 
			if dict_patient[key].max_death != 'unknown':
				dict_patient[key].disease_free_survival = min(dict_patient[key].min_new_tumor,dict_patient[key].max_death)
			else:
				dict_patient[key].disease_free_survival = dict_patient[key].min_new_tumor
		else:
			dict_patient[key].disease_free_survival = dict_patient[key].max_death
	else:
		if dict_patient[key].min_new_tumor != 'unknown': 
			if dict_patient[key].max_followup != 'unknown':
				dict_patient[key].disease_free_survival = min(dict_patient[key].min_new_tumor,dict_patient[key].max_followup)
			else:
				dict_patient[key].disease_free_survival = dict_patient[key].min_new_tumor
		else:
			dict_patient[key].disease_free_survival = dict_patient[key].max_followup

	if temp == 'Dead':
		dict_patient[key].censored = 1
	else:
		dict_patient[key].censored = 0


    

# write csv

headings = ['patient_id',
			'days_to_last_follow_up',
			'days_to_death',
			'days_to_new_tumor_event_after_initial_treatment',
			'vital_status',
			'days_to_last_known_alive',
			'tp53_mutation',
			'tp53_expression',
			'max_followup',
			'min_new_tumor',
			'valid_varient_classification',
			'max_death',
			'last_vital_status',
			'overall_survival',
			'disease_free_survival',
			'censored']
			

rows = []

def write_row(row,i,array):
	str1 = ''
	for index in range(len(array)):
		str1 += str(array[index])
		if index != len(array)-1:
			str1 += '|'
	row[i] = str1

for key in dict_patient:
	row = [None for x in range(20)]
	row[0] = key
	write_row(row,1,dict_patient[key].days_to_last_followup)
	write_row(row,2,dict_patient[key].days_to_death)
	write_row(row,3,dict_patient[key].days_to_new_tumor_event_after_initial_treatment)
	write_row(row,4,dict_patient[key].vital_status)
	write_row(row,5,dict_patient[key].days_to_last_known_alive)
	write_row(row,6,dict_patient[key].tp53_mutation)
	row[7] = dict_patient[key].tp53_expression
	row[8] = dict_patient[key].max_followup
	row[9] = dict_patient[key].min_new_tumor
	row[10] = dict_patient[key].valid_varient_classification
	row[11] = dict_patient[key].max_death
	row[12] = dict_patient[key].last_vital_status
	row[13] = dict_patient[key].overall_survival
	row[14] = dict_patient[key].disease_free_survival
	row[15] = dict_patient[key].censored
	rows.append(row)
	pprint(row)

with open('data.csv', 'w') as f:
	f_csv = csv.writer(f)
	f_csv.writerow(headings)
	f_csv.writerows(rows)
	




