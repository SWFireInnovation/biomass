__author__ = "Greg Cohn"
__email__ = "greg.cohn@nau.edu"
__version__ = "0.0"

import numpy as np
import yaml
import pandas as pd
import pint_pandas
import pint


class UnitsMngr:
    """
    units:
        in:
            DBH: cm
            HT: m
        out:
            weight: kg
        equ:
            DBH: in
            HT: ft
            weight: lbs
    """
    def __init__(self, df, units):
        self.df = df
        self.units = units

    def has_unit(self, series):
        return hasattr(series.dtype, 'unit') and series.dtype.units is not None

    def set_unit(self, col, unit):
        self.df[col] = self.df[col].astype(f'pint[{unit}]')

    def check_df_units(self, workflow):

        for col, unit in self.units[workflow].items():
            if not self.has_unit(self.df[col]):
                self.set_unit(col, unit)


class LoadData(UnitsMngr):
    def __init__(self, data_path_yaml='../data_info.yaml'):
        self.df = pd.DataFrame()

        self.units = {}
        self.data_path = ''
        self.set_data_info(data_path_yaml)

        self.required_columns = {'spp': {'alternate':['species'],
                                         'exclude':[]}}

        UnitsMngr.__init__(self, self.df, self.units)


    @staticmethod
    def load_yaml(filename):
        with open(filename, "r") as stream:
            return yaml.safe_load(stream)

    def set_data_info(self, data_loc_yaml='../data_info.yaml'):

        config =  self.load_yaml(data_loc_yaml) if '.yaml' in data_loc_yaml or '.yml' in data_loc_yaml else data_loc_yaml
        self.data_path = config['data_path']
        self.units =config['units']

    def load_data(self, filename):

        if filename.endswith(".csv"):
            self.df = pd.read_csv(filename)
        elif filename.endswith(".xlsx"):
            self.df = pd.read_excel(filename)

        self.check_req_columns(self.required_columns)

    def check_req_columns(self, required_columns):
        has_cols, rename_cols = self.find_req_column(required_columns)
        if not all(has_cols):
            missing = np.array(list(required_columns))[~has_cols]
            raise ValueError(f'Missing required columns.\nCannot find column {missing}')
        if rename_cols:
            self.df.rename(columns=rename_cols, inplace=True)

    def find_req_column(self, required_columns):
        # transform required columns into a list with all possible alternative column names
        dfcols_lower = np.array([c.lower() for c in self.df.columns])
        dfcols_lower_map = {c.lower(): c for c in self.df.columns}

        has_cols = []
        rename_cols = {}
        for req_col in required_columns:
            req_col_lwr = req_col.lower()
            # is required column in columns?
            if req_col in self.df.columns:
                has_cols.append(True)
                continue
            # remove case sensitivity and look for the column again
            elif req_col_lwr in dfcols_lower:
                orig_name = dfcols_lower_map[req_col_lwr]
                rename_cols[orig_name] = req_col
                has_cols.append(True)
                continue

            # check if the required column or a similar alternative is contained within a larger column name,
            # for example 'dbh' is within the column name 'dbh/drc' or check for 'species' instead of 'spp'
            alt = required_columns[req_col].get('alternate', [])
            candidates = [req_col_lwr] + ([alt] if isinstance(alt, str) else alt)
            include = np.zeros(len(dfcols_lower), dtype=bool)
            for req in candidates:
                include |= [req.lower() in col for col in dfcols_lower]

            exclude =  np.zeros(len(dfcols_lower), dtype=bool)
            for ex in required_columns[req_col].get('exclude', []):
                exclude |= [ex.lower() in col for col in dfcols_lower]

            orig_name_lwr = dfcols_lower[include & ~exclude]
            if orig_name_lwr.size > 0:
                orig_name = dfcols_lower_map[orig_name_lwr[0]]
                rename_cols[orig_name] = req_col
                has_cols.append(True)
            else:
                has_cols.append(False)


        return np.array(has_cols), rename_cols
