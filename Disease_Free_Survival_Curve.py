import csv
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from pprint import pprint
from matplotlib import pyplot as plt

T = []
E = []
T_valid = []
E_valid = []
T_invalid = []
E_invalid = []


with open('data.csv','r') as f:
	f_csv = csv.DictReader(f)
	for row in f_csv:
		if row['disease_free_survival'] != 'unknown':
			T.append(int(row['disease_free_survival']))
			if row['min_new_tumor'] == 'unknown':
				E.append(0)
			else:
				E.append(1)
			if row['valid_varient_classification'] == 'Yes':
				T_valid.append(T[-1])
				E_valid.append(E[-1])
			else:
				T_invalid.append(T[-1])
				E_invalid.append(E[-1])
'''
pprint(T)
pprint(E)
'''
print('T_valid\tE_valid')
for i in range(len(T_valid)):
	print(T_valid[i],'\t',E_valid[i])

results = logrank_test(T_valid, T_invalid, E_valid, E_invalid, alpha=0.95)
results.print_summary()

kmf = KaplanMeierFitter()

kmf.fit(T_valid, event_observed=E_valid, label='TP53 mutant')

kmf.survival_function_
kmf.median_
ax = kmf.plot(show_censors=True)




kmf.fit(T_invalid, event_observed=E_invalid, label='TP53 Normal')
kmf.plot(ax=ax,show_censors=True)
plt.title('Kaplan-Meier_Disease_Free_Survival_Curve')
plt.show()

