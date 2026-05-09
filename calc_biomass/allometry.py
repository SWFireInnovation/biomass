__author__ = "Greg Cohn"
__email__ = "greg.cohn@nau.edu"
__version__ = "0.0"

import numpy as np
import pandas as pd
from pint import UnitRegistry

ureg = UnitRegistry()

def _vect_by_spp(trees, spp_coef, equ):
    """
    Vectorize by species.

    For each species in the dictionary spp_coef, it selects that species from the pd.DataFrame trees and applies the
    equations with the coefficients supplied by spp_coef. If the species in spp_coef do not match any species in trees,
    it will simply return all 0's.

    :param trees: pd.DataFrame containing tree data with the columns necessary for the equation provided and 'spp'.
    :param spp_coef: dict where the keys are species and the values are coefficients.
    :param equ: method object that takes trees and a list of coefficients as inputs.
    :return: a list containing the output of equ.
    """
    out = np.zeros(len(trees))

    for spp, coef in spp_coef.items():
        this_spp = trees['spp_eq'] == spp
        out[this_spp] = equ(trees[this_spp], coef)

    return out

def eq_VOLUME_cylinder(trees, units):
    """
    Calculate the volume of a cylinder.

    Outputs are in m3

    :param trees: a pd.DataFrame from a class's `self.trees'
    :return: a pd.Series of tree volume in liters calculated as a cylinder for the slice `neg_vol`.
    """

    ht = trees['HT'].to_numpy() * eval(f'ureg.{units["HT"]}')
    ht_m = ht.to('m')

    dbh = trees['DBH'].to_numpy() * eval(f'ureg.{units["DBH"]}')
    dbh_m = dbh.to('m')

    a_m2 = np.pi * (dbh_m/2)**2
    v_m3 = a_m2 * ht_m

    return v_m3.magnitude

class OakWoodland_Chojnacky:
    """
    This class contains biomass equations developed based on sampling of a number of oak species common in AZ and NM.
    These equations focus on estimating total wood volume, and deriving estimates of foliage and small branchwood from
    total wood volume.

    Diameter at Root Collar (DRC) required for all equations.

    Input pd.DataFrame must be unit aware using `pint_pandas`.

    Chojnakcy, D.C. 1992 . Estimating volume and biomass for dryland oak speices. In: Ffolliott P.F., Gottfried, G.J.,
    Bennett, D.A., Hernandez, C.V.-M., Ortega-Rubio, A., and R.H. Hamre, technical coordinators. Ecology and management
    of oak and associated woodlands: perspectives in the southwestern United States and Northern Mexico. Proceedings
    April 7-30, 1992; Sierra Vista, AZ. USDA Forest Service, Rocky Mountain Forest and Range Experiment Station
    General Technical Report RM-218:155-161.

    Species used in the development of this equation include: Quercus arizonica; Q. emoryi; Q. gambelii;
    Q. hypoleucoides; Q. oblongifolia, Q. ilex.

    """
    units = {'equ':{'DBH': 'cm',
                    'HT':'m'},
            'out':{'weight': 'kg'}
            }
    required_columns = {'DBH': {'alternate': ['drc'],
                                'exclude': ['circumfrence', 'circumference']},
                         'HT':{'alternate':['height'],
                               'exclude':['canopy', 'base', 'lowest']}
                        }
    def __init__(self, trees):
        self.trees = trees.fillna(0)

    def correct_negative_vol(self, vol_bd):
        """
        For very small trees (<7.5 cm DRC, <4 m HT) :meth:`eq_VOLUME_branchdiam` yields negative numbers. In these
        instances, a simple volume of a cylinder is applied using height and dbh.

        :param vol_bd: list of volume in liters
        :return: list of volume in liters
        """
        neg_vol = vol_bd < 0
        # 5th percentile of DRC for species in original study was 8 -12 cm.
        # 7.5 was taken from visual assessment of where the equ could reasonably produce a negative number
        neg_dbh = self.trees['DBH'] < 7.5
        # also taken from visual assessment
        neg_ht = self.trees['HT'] < 4

        neg_eq = neg_vol & neg_dbh & neg_ht

        # output is in m3, convert to liters (dm3)
        vol_bd[neg_eq] = eq_VOLUME_cylinder(self.trees[neg_eq], self.units['equ'])*1000

        # if there is a negative volume with a large dbh or a large height, the data probably wasn't entered correctly
        error = neg_vol & ~(neg_dbh & neg_ht)
        if error.any():
            raise ValueError(f'Negative tree volume!\n{self.trees[error]}')
        else:
            return vol_bd

    def eq_VOLUME_branchdiam(self):
        """
        Equation to calculate total tree volume of branches greater than 3.8 cm.

        This equation uses DRC to estimate the total tree volume of tree bole and branches >3.8cm. The training data
        was derived from the sum all measured branch and stem segments. Outputs are in liters (decimeters**3). Only the
        combined equation for both single and multi-stem trees is provided below.

        Species used in the development of this equation include:
        Quercus arizonica; Q. emoryi; Q. gambelii; Q. hypoleucoides; Q. oblongifolia.

        :return: list of tree volume in decimeters**3 (liters).
        """

        # Uses combined multiple stem and single stem equation
        # breakpoint between equation forms (metric)
        X0 = 6.9850
        # equation parameters (metric)
        B0, B1, B2 = -1.5107, 66.462, 1.3694

        # DRC**2 * HT / 2000 metric independent variable
        X = (self.trees['HT'] * self.trees['DBH']**2)/2000

        equ1 = X < X0
        equ2 = X > X0

        # empty list of volume for branchdiameter
        vol_bd = np.zeros(len(X))

        vol_bd[equ1] = B0 + B1*X[equ1] + B2*X[equ1]**2
        vol_bd[equ2] = B0 + B1*X[equ2] + B2*(3*X0**2 - 2*X0**3/X[equ2])

        return self.correct_negative_vol(vol_bd)

    @staticmethod
    def _eq_WEIGHT_branchdiam(volbd, coef):
        sg, Dh20 = coef
        return volbd['vol_bd'] * sg * Dh20

    def eq_WEIGHT_branchdiam(self, vol_bd):
        """
        Equation to calculate stem weight from volume.

        Takes volume of stem and branchwood and converts to weight by multiplying by specific gravity. Volume is
        assumed to be in decimeters**3 (liters), so a constant of 1 kg per liter is used.

        Original paper used presented specific gravity of Emory Oak and Gambel oak. Most other species must use a
        genus level average from Miles and Smith 2009 (RN-NRS-38).

        :param vol_bd: list of stem volume in decimeters**3 (liters).
        :return: list of weight by tree.
        """
        # density of water kg/liter  (kg/decimeter**3)
        Dh20 = 1
        spp_sg = {# from Chojnacky 1992 + Maingi and Ffolliott 1992
                  # emory and gambel oak
                  'QUGA': [0.6340, Dh20],
                  'QUEM': [0.5670, Dh20],
                  # from Chojnacky, et. al. 2013
                  # silver leaf oak
                  'QUHY': [0.59, Dh20],
                  # Miles and Smith 2009 (RN-NRS-38)
                  # generic average of oak species
                  'QUERCUS': [0.59, Dh20],
                  }

        # slicing by species requires a species column!
        spp_vol_bd = pd.DataFrame({'spp_eq':self.trees['spp_eq'], 'vol_bd':vol_bd})
        return _vect_by_spp(spp_vol_bd, spp_sg, self._eq_WEIGHT_branchdiam)

    def _eq_wght_crwn(self, wght_bd, coefs):
        B0, B1, B2 = coefs
        return 10**(B0 + B1*np.log10(wght_bd) + B2*self.trees['HT'].to_numpy())

    def eq_WEIGHT_foliage(self, wght_bd):
        """
        Equation to calculate foliage weight from stem weight and height.

        :param wght_bd: list of stem weight in Kg.
        :return: list of foliar weight in Kg.
        """
        coef = [-0.6210, 0.8382, -0.0307]

        return self._eq_wght_crwn(wght_bd, coef)

    def eq_WEIGHT_branch(self, wght_bd):
        """
        Equation to calculate weight of branches < 5cm from stem weight and height.

        .. warning::
            :meth:`eq_WEIGHT_branchdiam` calculates weight of all stem sections down to 3.8 cm. The weight of branch
            tips calculated by this equation will dupliacte weights for branches 3.8-5cm in diameter. Chojnakcy felt
            this was likely to be <10% overlap.

        :param wght_bd: list of stem weight in Kg.
        :return: list of branch weight in Kg.
        """
        coef = [0.2264, 0.7752, -0.0132]

        return self._eq_wght_crwn(wght_bd, coef)

    def calc_tot_tree_weight(self):
        wght_bole =self.eq_WEIGHT_branchdiam(self.eq_VOLUME_branchdiam())
        wght_br = self.eq_WEIGHT_branch(wght_bole)
        wght_fol = self.eq_WEIGHT_foliage(wght_bole)

        return wght_bole + wght_br + wght_fol

    def calc_foliar_weight(self):
        wght_bole = self.eq_WEIGHT_branchdiam(self.eq_VOLUME_branchdiam())

        return self.eq_WEIGHT_foliage(wght_bole)

    def calc_avl_canfuel(self):
        """
        Calculate the canopy fuel load as the proportion of the crown available to the advancing front of a crown fire;
        defined as foliage + 1/2 1hr fuels.

        .. Warning::
            This does not represent a true canopy fuel load because the branch does not represent the 1 hr fuel class.
            The branch equation calculates weight of all branches up to 5cm.

        :return: list of total available canopy fuel in Kg.
        """
        wght_bole = self.eq_WEIGHT_branchdiam(self.eq_VOLUME_branchdiam())
        wght_br = self.eq_WEIGHT_branch(wght_bole)
        wght_fol = self.eq_WEIGHT_foliage(wght_bole)

        return wght_br/2 + wght_fol

class PJ_Grier:
    """
    These equations were developed to estimate biomass of pinyon and juniper on the Mogollon rim in AZ outside of
    Flagstaff.Sample sizes of destructive sampling were small (<=15).

    Diameter at Root Collar (DRC) is required for all equations.

    Input pd.DataFrame must be unit aware using `pint_pandas`.

    To get accurate estimates of the canopy fuel available to the moving front of a crown fire, an additional
    coefficient is used to adjust the branchwood weight to represent 1/2 of the 1 hour fuels following the style of
    FuelCalc.

    Grier, C.C., Elliott, K.J., McCullough, D.G. 1992. Biomass distribution and productivity of Pinus edulis-Juniperus
    monosperma woodlands of north-central Arizona. Forest Ecology and Management, 50:331-350.
    """
    units = {'equ':{'DBH': 'cm'},
            'out':{'weight': 'kg'}
            }
    required_columns = {'DBH': {'alternate': ['drc'],
                                'exclude': ['circumfrence', 'circumference']}
                             }

    def __init__(self, trees):
        self.trees = trees

    @staticmethod
    def eq_weight(DRC, coef):
        """
        Logarithmic equation used to calculate species component weights.

        .. Note::
            Adds an additional parameter frac. Frac is an adjustment factor used by FuelCalc to represent the
            proportion of branchwood that is equivalent to 0.5*1hr, since this equation provides weight of 1 and 10hr's.

        Returns weight in Kg.

        :param DRC: pd.DataFrame with a column labeled 'DBH' containing DRC's.
        :param coef: list of a and b coefficients for equation.
        :return: list of weight in Kg.
        """
        a, b, frac = coef
        return frac * 10**(a + b * np.log10(DRC['DBH']))

    def calc_tot_tree_weight(self):
        """
        Calculate total tree weight in Kg.

        .. Note::
            Additional coefficient multiply output by 1.

        :return: list of total tree weight in Kg.
        """

        coef = {"PIED": [-1.468, 2.582, 1],
                "JUMO": [-1.157, 2., 1]
                }

        return _vect_by_spp(self.trees, coef, self.eq_weight)
        #self._calc_spp_weight(coef)

    def calc_foliar_weight(self):
        """
        Calculate foliar weight in Kg.

        .. Note::
            Additional coefficient multiply output by 1.

        :return: list of foliar weight in Kg.
        """
        coef = {"PIED": [-0.946, 1.565, 1],
                "JUMO": [-1.737, 1.382, 1]
                }

        return _vect_by_spp(self.trees, coef, self.eq_weight)
        #self._calc_spp_weight(coef))

    def calc_1hr_weight(self):
        """
        Calculate the estimated weight of 1/2 of branches in 1hr size class.

        The original equation calculates the weight of all branches up to 2.5 cm. Adjustment fractions are taken from
        FuelCalc to represent the proportion of the weight that is 1/2 of the 1hr size class.

        :return: list of weight in Kg.
        """
        branch_coef = {"PIED": [-1.613, 2.088, 0.33],
                       "JUMO": [-1.476, 1.787, 0.25]
                       }
        return _vect_by_spp(self.trees, branch_coef, self.eq_weight)

    def calc_avl_canfuel(self):
        """
        Calculate the weight of canopy fuel available to the flaming front of a crown fire in Kg (foliage + 0.5*1hr fuel).

        .. Warning::
            Branch biomass includes all branches up to 2.5 cm (10 hour fuels). This equation over-estimates canopy fuel.
            To adjust, the branchwood is reduced in the style of FuelCalc by adding a species specific
            adjustment to the coefficients. See :meth:`calc_1hr_weight`.

        :return: list of available canopy fuel in Kg.
        """

        wght_fol = self.calc_foliar_weight()
        wght_br = self.calc_1hr_weight()

        return wght_br + wght_fol

class BCtimber_Standish:
    """
    This class contains biomass equations developed for major species found in British Columbia. Separate equations are
    provided for wood, bark, 3 branch size classes, and foliage. All coefficients are additive, so that the sum of all
    coefficients will calculate total tree biomass.

    Requires DBH, height, and volume as inputs. Volume equations are included in this class and are taken from the BC
    Ministry of Forestry volume tables.

    Input pd.DataFrame must be unit aware using `pint_pandas`.

    Standish, J.T., Manning, G.H., and Demaerschalk, J.P. 1985. Development of biomass equations for British Columbia
    tree species. Canadian Forestry Service, Pacific Forest Research Centre, Information Report BC-X-264 (Vancouver, BC)

    Susan Watts eds. 1983. Forestry Handbook for British Columbia, 4th Edition. University of British Columbia. Forest
    Club.
    """
    units = {'equ':{
                    'HT': 'meters',
                    'DBH': 'meters'
                    },
            'out':{'weight': 'kg'}
            }
    required_columns = {'DBH': {'alternate': ['drc'],
                                'exclude': ['circumfrence', 'circumference']},
                         'HT': {'alternate': ['height'],
                                'exclude': ['canopy', 'base', 'lowest']}}

    def __init__(self, trees):
        self.trees = trees

        self.vol = None

    def correct_vol(self, vol):
        """
        Correct for disproportionate volumes that will lead to negative biomass.

        It was found that dead snags ("candle sticks") can have negative biomass. This occurs because the multiple
        regression model can combine a term for volume with a term for something like diameter squared. If the tree is
        very short, but very thick (i.e. broken topped snag), the volume equation underestimates and cannot balance the
        diameter term. This is assumed to be a snag, where a cylinder is a more appropriate approximation of volume.

        :param vol: a pd.Series of tree volumes calculated by :meth:`BCtimber_Standish.calc_vol()`.
        :return:
        """

        stubby = (self.trees['HT'] / self.trees['DBH']) < 8

        cyl = eq_VOLUME_cylinder(self.trees, self.units['equ'])
        undervol = (vol/cyl) < 0.365

        short = self.trees['HT'] < 4

        broken_snag = stubby & short & undervol

        vol[broken_snag] = cyl[broken_snag]

        return vol

    def equ_weight(self, tr, coef):
        """
        Equation for tree biomass. For any given species, a maximum of 3 terms are used. Coefficients are 0 for any
        unused terms.

        :param tr: pd.DataFrame containing columns DBH, height, and volume
        :param coef: list of 8 coefficients for this species.
        :return: list of component weight in Kg.
        """

        b0, b_v, b_d2, b_hv, b_d2hv, b_dhv, bd2v, d = coef

        wt = b0 + b_v*tr['vol'] + \
             b_d2*tr['DBH']**2 + \
             b_hv*tr['HT']*tr['vol'] +\
             b_d2hv*tr['vol']*tr['HT']*tr['DBH']**2 + \
             b_dhv*tr['vol']*tr['HT']*tr['DBH'] + \
             bd2v*tr['vol']*tr['DBH']**2 + \
             d*tr['DBH']

        return wt

    def equ_vol(self, tr, coef):
        """
        Equation for bole volume from British Columbia Ministry of Forestry volume tables.

        .. Note::
            DBH for all biomass equations are in m.

        Forest Inventory Division. 1976. Whole stem cubic metre volume tables: centimetre diameter class merchantable
        volume factors. B.C. Forest Service.In: Susan Watts eds. 1983. Forestry Handbook for British Columbia: 4th ed.
        Forestry Undergraduate Society, Faculty of Forestry. University of British Columbia.

        :param tr: pd.DataFrame containing columns DBH and HT.
        :param coef: A list of 3 coefficients.
        :return: list of tree volume in m3.
        """

        b0, b1, b2 = coef

        dbh_cm = tr['DBH'] * 100

        v = b0 + b1*np.log10(dbh_cm) + b2*np.log10(tr['HT'])

        return 10**v

    def calc_vol(self):
        """
        Calculate tree volume using Ministry of Forestry volume tables from BC, Canada.

        Susan Watts eds. 1983. Forestry handbook for british columbia, 4th Edition.

        :return:
        """
        coef = {'POTR5':[-4.419728, 1.894760, 1.053730],
                'PICO':[- 4.349504, 1.822760, 1.108120],
                'ABLA':[- 4.291919, 1.872930, 0.998274], # all "balsam" species the same
                'PIEN':[- 4.294193, 1.858590, 1.007790], # all spruce species the same
                'PIPO':[-4.482485, 1.954430, 1.01677],
                'PSME':[-4.383102, 1.742940, 1.15641],
                'POTR15':[-4.648431, 1.735180, 1.356010], # %all cottonwood species the same
                'ABGR':[-4.291919, 1.872930, 0.998274], # %all "balsam" species the same.
                'LAOC':[-4.350486, 1.723600, 1.135270],
                'TSHE':[-4.394633, 1.942900, 0.990275],
                'THPL':[-4.178431, 1.759950, 1.019080],
                'BEPA':[-4.443142, 1.909560, 1.052050],
                'PIMO3':[-4.300522, 1.857800, 1.022250],
                'TSME':[-4.394633, 1.942900, 0.990275] #; %TSHE substituted for TSME
                }

        vol = _vect_by_spp(self.trees, coef, self.equ_vol)
        return self.correct_vol(vol)


    def set_vol(self):
        """
        Calculate tree bole volume and assign to class instance.
        """
        self.vol = self.calc_vol()

    def calc_tot_tree_weight(self):
        coef = {'POTR5':[-10, 311, 2205, 0, 0, 0, 0, 0],
                'PICO':[-12.9, 139.5, 3619.7, 0, 0, 0, 0, 0,],
                'ABLA':[47.8, - 106, 8040.9, 0, 0, 0, 0, -1093.6],
                'PIEN':[17.6, 415.9, 611.5, 0, 0, 0, 0, 0],
                'PIPO':[23.9, 733.1, -851.9, 0, 0, 0, 0, 0],
                'PSME':[-28.5, -242.4, 5800.3, 0, 4.1, 0, 0, 0],
                'POTR15':[2, 460.9, 695.6, 0, 0, 0, 0, 0],# %Populus trichocarpa substituted for POAN
                'ABGR':[-1.7, 270.9, 1658.3, 0, 0, 0, 0, 0],
                'LAOC':[-24, 1395, 0, -23.8, 0, 0, 0, 0],
                'TSHE':[-17, 756.3, 0, 0, 0, 0.7, -518.7, 0],
                'THPL':[13.3, 320.3, 539, 0, 0, 0, 0, 0],
                'BEPA':[-6.5, 755.3, 348.7, 0, 0, 0, 0, 0],
                'PIMO3':[-10.3, 334.2, 1628.2, 0, 0, 0, 0, 0],
                'TSME':[11.5, 642.2, 0, 0, -5.9, 464.1, 0, 0]
                }

        if self.vol is None:
            self.set_vol()

        tree = self.trees.loc[:, ['DBH', 'HT', 'spp_eq']]
        tree['vol'] = self.vol

        return _vect_by_spp(tree, coef, self.equ_weight)

    def calc_foliar_weight(self):
        coef = {'POTR5':[0.5, 0.2, 89.2, 0, 0, 0, 0, 0],
                'PICO':[-1.3, -41.5, 569.2, 0, 0, 0, 0, 0],
                'ABLA':[15.2, -175.1, 2614.1, 0, 0, 0, 0, -332.6],
                'PIEN':[7.2, -18.4, 348.5, 0, 0, 0, 0, 0],
                'PIPO':[0.2, -18.9, 392.9, 0, 0, 0, 0, 0],
                'PSME':[-3.4, -79.5, 860.0, 0, -0.6, 0, 0, 0],
                'POTR15':[0.2, -0.4, 156.6, 0, 0, 0, 0, 0],#        %Populus trichocarpa substituted for POAN
                'ABGR':[-2.9, -37, 757.4, 0, 0, 0, 0, 0],
                'LAOC':[-0.3, 55.9, 0, -1.4, 0, 0, 0, 0],
                'TSHE':[2.7, 43.3, 0, 0, 0, 2.6, -220.9, 0],
                'THPL':[-4, -42.4, 554.6, 0, 0, 0, 0, 0],
                'BEPA':[1, 6.9, 119.3, 0, 0, 0, 0, 0],
                'PIMO3':[-1, 36.9, 575.3, 0, 0, 0, 0, 0],
                'TSME':[1.8, 57.1, 0, 0, -19.3, 914.5, 0, 0]
                }

        if self.vol is None:
            self.set_vol()

        tree = self.trees.loc[:, ['DBH', 'HT', 'spp_eq']]
        tree['vol'] = self.vol

        return _vect_by_spp(tree, coef, self.equ_weight)

    def calc_1hr_weight(self):
        coef = {'POTR5':[-1.2, -29, 271.1, 0, 0, 0, 0, 0],
                'PICO':[-0.7, -9.3, 153.7, 0, 0, 0, 0, 0],
                'ABLA':[5.4, -99.5, 1457.5, 0, 0, 0, 0, -165.4],
                'PIEN':[0.2, -23.6, 317.3, 0, 0, 0, 0, 0],
                'PIPO':[0.3, -0.5, 5.8, 0, 0, 0, 0, 0],
                'PSME':[-0.3, -99.1, 628.9, 0, 4.3, 0, 0, 0],
                'POTR15':[-0.1, -7.6, 80.3, 0, 0, 0, 0, 0],# %Populus trichocarpa substituted for POAN
                'ABGR':[0.2, -10.3, 193.5, 0, 0, 0, 0, 0],
                'LAOC':[-1.3, 76.3, 0, -1.9, 0, 0, 0, 0],
                'TSHE':[-0.1, 33.1, 0, 0, 0, 1.4, -137.5, 0],
                'THPL':[0.7, -12.2, 130.1, 0, 0, 0, 0, 0],
                'BEPA':[0.4, 23.4, -3, 0, 0, 0, 0, 0],
                'PIMO3':[-0.5, 6.7, 61.7, 0, 0, 0, 0, 0],
                'TSME':[0.6, -4.3, 0, 0, -7.7, 561.7, 0, 0]
                }

        if self.vol is None:
            self.set_vol()

        tree = self.trees.loc[:, ['DBH', 'HT', 'spp_eq']]
        tree['vol'] = self.vol

        return _vect_by_spp(tree, coef, self.equ_weight)

    def calc_avl_canfuel(self):
        if self.vol is None:
            self.set_vol()

        wght_br = self.calc_1hr_weight()
        wght_fol = self.calc_foliar_weight()

        return wght_br/2 + wght_fol


class Mesquite_McClaren:
    """
    This class contains biomass equations developed for velvet mesquite trees (Prosopis velutina) in southern AZ,
    south of Tucson. The sample size is small, but coefficient b=2.19 in the equation ln(y) = a+b*lnX where X is total
    biomass is within the range of similar findings for other mesquite species where b ranged from 2.1-2.37
    (Alvarez et al., 2011; Northup et al. 2005; Padron and Navarro, 2004). Navar et al. 2019 also found that this
    equation was within the lower range of data of the 510 samples in Navar et al. 2019.

    It is also important to note that the "fine stem" category of branchwood is 0-1cm. The standard "available canopy
    fuel" is foliage + 0.5* 1hr's. 1hr's are 0-6 mm, so the fine stem category will over-estimate canopy fuel.

    Diameter at Root Collar (DRC) required for all equations.

    Input pd.DataFrame must be unit aware using `pint_pandas`.

    McClaran, M.P., McMurtry, C.R., Archer, S.R.. 2013. A tool for estimating impacts of woody encroachment in arid
    grasslands: Allometric equations for biomass, carbon and nitrogen content in Posopis veluntina. Journal of Arid
    Environnments 88(2013): 39-12
    """
    units = {'equ':{'DBH': 'cm'},
            'out':{'weight': 'kg'}
            }
    required_columns = {'DBH': {'alternate': ['drc'],
                                'exclude': ['circumfrence', 'circumference']}
                        }

    def __init__(self, trees):
        self.trees = trees

    @staticmethod
    def eq_weight(DRC, coef):
        """
        Logarithmic equation (natural) used to calculate species component weights.

        Returns weight in Kg.

        :param DRC: pd.DataFrame with a column labeled 'DBH' containing DRC's.
        :param coef: list of a and b coefficients for equation.
        :return: list of weight in Kg.
        """
        a, b, CF = coef
        return np.e ** (CF * (a + b * np.log(DRC['DBH'])))

    def calc_tot_tree_weight(self):
        """
        Calculate total tree weight in Kg.

        :return: list of total tree weight in Kg.
        """

        coef = {"PRVE": [-3.02, 2.19, 1.06]}

        return _vect_by_spp(self.trees, coef, self.eq_weight)

    def calc_foliar_weight(self):
        """
        Calculate foliar weight in Kg.

        :return: list of foliar weight in Kg.
        """
        coef = {"PRVE": [-4.88, 1.67, 1.02]}

        return _vect_by_spp(self.trees, coef, self.eq_weight)
        # self._calc_spp_weight(coef))

    def calc_1hr_weight(self):
        """
        Calculate the estimated weight of branches in 1hr size class.

        It is important to note that the "fine stem" category of branchwood is 0-1cm. 1 hr fuel classes are 0-6mm, so
        the fine stem category will over-estimate canopy fuel.

        :return: list of weight in Kg.
        """
        coef = {"PRVE": [-3.15, 1.52, 1.01]}

        return _vect_by_spp(self.trees, coef, self.eq_weight)

    def calc_avl_canfuel(self):
        """
        Calculate the weight of canopy fuel available to the flaming front of a crown fire in Kg (foliage + 0.5*1hr fuel).

        .. Warning::
            Branch biomass includes all branches up to 1 cm. This equation over-estimates canopy fuel.

        :return: list of available canopy fuel in Kg.
        """

        wght_fol = self.calc_foliar_weight()
        wght_br = self.calc_1hr_weight()

        return wght_br/2 + wght_fol
