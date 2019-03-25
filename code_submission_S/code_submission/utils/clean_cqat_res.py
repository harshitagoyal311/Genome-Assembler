import pandas as pd
import sys


if __name__ == '__main__':
	#Reading the file  --- Temporarily not required
	df = pd.read_csv(sys.argv[1], sep='\t')
	cols = list(df.columns)
	req_cols = [cols[0], cols[5]]
	new_df = df[req_cols]
	new_df.columns = ['ID', 'Quality_Score']
	new_df.to_csv(str(sys.argv[1].split('.')[0] +'_cleaned.csv'), index=None )
