__author__ = "Greg Cohn"
__email__ = "greg.cohn@nau.edu"
__version__ = "0.0"

import numpy as np
import yaml
import pandas as pd
import pint_pandas
import pint
from calc_biomass import allometry

# unify the units registry so that we can check units
ureg = pint.UnitRegistry()
pint_pandas.PintType.ureg = ureg


class UnitsMngr:
    """
    units:
        input:
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

    def has_units(self, series):
        return hasattr(series.dtype, 'units') and series.dtype.units is not None

    def has_same_unit(self, series, target_unit):
        return series.dtype.units == ureg.Unit(target_unit)

    def assign_unit(self, col, unit):
        return self.df[col].astype(f'pint[{unit}]')

    def convert_unit(self, col, unit):
        return self.df[col].pint.to(unit)

    def get_df_units(self, workflow):

        df = pd.DataFrame(index=self.df.index)
        for col, unit in self.units[workflow].items():
            series = self.df[col]
            if not self.has_units(series):
                df[col] = self.assign_unit(col, unit)
            elif not self.has_same_unit(series, unit):
                df[col] = self.convert_unit(col, unit)
            elif self.has_units(series) and self.has_same_unit(series, unit):
                df[col] = series

        return df

    def check_df_units(self, workflow):

        cols = list(self.units[workflow])
        self.df[cols] = self.get_df_units(workflow)



class LoadData(UnitsMngr):
    def __init__(self, data_path_yaml='../data_info.yaml'):
        self.df = pd.DataFrame()

        # load the following from yaml
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

        self.df = self.check_req_columns(self.df, self.required_columns)
        self.check_df_units('input')

    def check_req_columns(self, df, required_columns):
        has_cols, rename_cols = self.find_req_column(df, required_columns)
        if not all(has_cols):
            missing = np.array(list(required_columns))[~has_cols]
            raise ValueError(f'Missing required columns.\nCannot find column {missing}')
        if rename_cols:
            df.rename(columns=rename_cols, inplace=True)

        return df

    def find_req_column(self,df, required_columns):
        # transform required columns into a list with all possible alternative column names
        dfcols_lower = np.array([c.lower() for c in df.columns])
        dfcols_lower_map = {c.lower(): c for c in df.columns}

        has_cols = []
        rename_cols = {}
        for req_col in required_columns:
            req_col_lwr = req_col.lower()
            # is required column in columns?
            if req_col in df.columns:
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


class CalcBiomass:
    def __init__(self, loaded_data, species_eq='../species_equ.yaml'):

        self.ldata = loaded_data

        self.species_eq = species_eq

    def get_spp_equ(self):
        return LoadData.load_yaml(self.species_eq)

    def add_equ_column(self):
        df = self.ldata.df
        spp_eq = self.get_spp_equ()

        for spp, eq in spp_eq.items():
            df.loc[df['spp']==spp, 'spp_eq'] = eq['spp_eq']
            df.loc[df['spp']==spp, 'equ'] = eq['equ']

    def get_clean_df(self, df, required_columns, units):
        # make sure there are the required named columns
        # find correct columns and adjust names if necesary
        df_req_col = LoadData().check_req_columns(df, required_columns)
        #convert to correct units and only save columns with defined units
        df_correct_units = UnitsMngr(df_req_col, units).get_df_units('equ')

        # remove unit data types and drop unit header
        # units cannot be used in complex equations
        df_clean = df_correct_units.pint.dequantify().droplevel(level=1, axis=1)
        # add species information
        df_clean[['spp_eq', 'equ']] = df[['spp_eq', 'equ']]

        return df_clean

    def calc_equ(self, data, equ, component):

        cls = getattr(allometry, equ)
        # check for the required columns and units of the equ
        df = self.get_clean_df(data, required_columns=cls.required_columns, units=cls.units)

        # run the equation
        cls_inst = cls(df)

        if component == 'total':
            wt = cls_inst.calc_tot_tree_weight()
        elif component == 'foliage' or component == 'foliar':
            wt = cls_inst.calc_foliar_weight()
        elif 'fuel' in component:
            wt = cls_inst.calc_avl_canfuel()

        # assign equ units to the output data (UnitsMngr needs a pd.DataFrame)
        wt = pd.DataFrame(data=wt, columns=['weight'])
        wt['weight'] = UnitsMngr(wt, cls_inst.units).assign_unit('weight', cls_inst.units['out']['weight'])

        # convert to desired output units
        mng_unit = UnitsMngr(wt, self.ldata.units)

        return mng_unit.get_df_units('out')

    def calc_biomass(self, component):

        equ_used = self.ldata.df.equ.unique().dropna()
        df = self.ldata.df

        wt_unit = self.ldata.units["out"]["weight"]
        df[f'{component}_wt'] = pd.Series(data=np.nan, index=df.index, dtype=f'pint[{wt_unit}]')
        for eq in equ_used:
            index_eq = df['equ'] == eq
            df.loc[index_eq, f'{component}_wt'] = self.calc_equ(df[index_eq], eq, component)['weight']


if __name__ == '__main__':
    data = LoadData()
    # wt = manage_data.CalcBiomass()

    data.load_data(data.data_path)
    data.df.rename(columns={'Final DBH/DRC': 'DBH', 'spp': 'Species', 'Species Code': 'spp'}, inplace=True)

    wt = CalcBiomass(data)
    wt.add_equ_column()
    wt.calc_biomass(component='foliar')

    wt.calc_biomass('total')