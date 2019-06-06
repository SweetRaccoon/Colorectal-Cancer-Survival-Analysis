import csv
from lifelines import KaplanMeierFitter
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
		if row['overall_survival'] != 'unknown':
			T.append(int(row['overall_survival']))
			if row['last_vital_status'] == 'Dead':
				E.append(1)
			else:
				E.append(0)
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

kmf = KaplanMeierFitter()

kmf.fit(T_valid, event_observed=E_valid, label='valid mutated')

kmf.survival_function_
kmf.median_
ax = kmf.plot()




kmf.fit(T_invalid, event_observed=E_invalid, label='invalid mutated')
kmf.plot(ax=ax)
plt.title('Kaplan-Meier Survival Curve')
plt.show()

