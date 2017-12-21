import numpy as np, astropy.io.fits as fits, matplotlib.pyplot as plt
from astropy.table import Table
from scipy.optimize import minimize
import emcee as mc

class Star:
    
    def __init__(self, ra, dec):
        """
        Parameters
        ----------
        ra : float
            Right ascension in decimal degrees
        dec : float
            Declination in decimal degrees
        """
        
        self.ra = ra
        self.dec = dec
        
        #Keeps track of how well covered this star is by our survey
        self.score = 0
        self.field_list = []
        
    def count_fields_star_in(self,test_field_list):
        """
        Counts how many fields this star is in.
        
        Parameter
        ---------
        test_field_list : list
            List of Field objects.
        
        Returns
        -------
        N_fields : int
            Number of fields this star is in.
        """
        
        N_fields = 0
        
        for field in test_field_list:
            
            if field.check_star_in_field(self):
                
                N_fields += 1
                self.field_list.append(field)
        
        return N_fields
    
    def score_star(self,overlap_bonus = 0.25):
        """
        If star is in one field, it gets 1 point. If its in more than one field, it gets 
        1 + overlap_bonus points. Must run count_fields_star_in first.
        
        Parameter
        ---------
        overlap_bonus : float
            bonus you want to give to stars that appear in multiple fields
        
        Returns
        -------
        score : float:
            This star's score!
        """
        
        if len(self.field_list) > 1:
            self.score = 1.0 + overlap_bonus
        elif len(self.field_list) == 1:
            self.score = 1
        else:
            self.score = 0
        
        return self.score
        
    def reset(self):
        """
        Resets field_list and score to empty for a new run
        """
        
        self.score = 0
        self.field_list = []
        
class Field:
    
    def __init__(self, ra, dec, size):
        """
        Parameters
        ----------
        ra : float
            Right ascension of center in decimal degrees
        dec : float
            Declination of center in decimal degrees
        size : float
            Size of square region on a side, arcminutes
        """
        
        self.ra = ra
        self.dec = dec
        self.size = size
    
    @classmethod
    def from_star_list(cls, star_list, size, random = True, i = 0):
        """
        Initialize field from list of stars. If random, will choose randomly. Otherwise
        it will just star_list[i]
        
        Parameters
        ----------
        star_list : list
            List of Star objects
        size : float
            Size of square region on a side, arcminutes
        random : bool
            If true, chooses randomly from star_list. Otherwise, chooses ith entry
        i : int
            If random is false, chooses star_list[i]
        
        """
        
        if random:
            star_init = np.random.choice(star_list)
        else:
            star_init = star_list[i]

        field = cls(star_init.ra,star_init.dec,size)

        return field
        
    def check_star_in_field(self,star):
        """
        Checks if a given star is in the field.
        
        Parameter
        ---------
        star : Star
            Star object
        
        Returns
        -------
        in_field : bool
            True if star is in the field, false if not
        """

        delta_dec_arcmin = np.abs(star.dec-self.dec) * 60.0

        delta_ra_arcmin = np.abs(star.ra - self.ra) * np.cos(star.dec * np.pi/180) * 60.0

        if (delta_ra_arcmin < self.size/2.0)&(delta_dec_arcmin < self.size/2.0):

            return True

        return False
    
    def count_stars_in_field(self,star_list):
        """
        Counts the number of stars in this field.
        
        Parameter
        ---------
        star_list : list
            List of Star objects.
        
        Returns
        -------
        N_in_field : int
            Number of stars in the field from star_list
        """
        
        N_in_field = 0
        
        for star in star_list:
            
            if self.check_star_in_field(star):
                
                N_in_field += 1
        
        return N_in_field
    
    def to_region_string(self, star_list = None):
        """
        Exports this field to a ds9 region
        
        Parameter
        ---------
        star_list : list
            list of Star objects. Counts the number of stars in the field, exports as text
        
        Returns
        -------
        region_string : str
            a ds9 compatible region string.
        """
        
        region_string = "box({0},{1},{2}',{2}',0)".format(self.ra,self.dec,self.size)
        
        if star_list is not None:
            
            region_string += ' # text={{}}'.format(self.count_stars_in_field(star_list))
            
        return region_string
    
def score_survey(star_list,field_list,overlap_bonus):
    """
    Given the stars you want and a set of fields, scores the stars, sums it up.
    
    Parameters
    ----------
    star_list : list
        list of Star objects
    field_list : list
        list of Field objects
    overlap_bonus
        
    Returns
    -------
    score : float
        Sum of the scores of all of the stars
    """
    
    score = 0.0
    
    for star in star_list:
        
        star.reset()
        
        star.count_fields_star_in(field_list)
        
        score += star.score_star(overlap_bonus = overlap_bonus)
        
    return score

def make_survey_and_score(theta,args):
    """
    Wrapper for score_survey. Makes a new field_list, then scores the stars in star_list.
    
    Parameter
    ---------
    theta : list
        Should be a list of length N_fieldx2, with first N_field entries being RAs, next 
        N_field entries being Decs.
    args : tuple
        First entry should be a list of star objects in your survey. Second entry
        should be the size of the field of view of the camera in arcminutes. Third entry 
        should be the overlap_bonus
        
    Returns
    -------
    survey_score : float
        Score of this particular survey
    """
    
    try:
        assert len(theta) % 2 == 0
    except AssertionError as e:
        raise AssertionError("Length of theta should be N_fields x 2!")
    
    n_fields = len(theta) // 2
    ras = theta[:n_fields]
    decs = theta[n_fields:]
    
    star_list = args[0]
    field_size = args[1]
    overlap_bonus = args[2]
    
    field_list = [Field(ra,dec,field_size) for ra,dec in zip(ras,decs)]
    
    return score_survey(star_list,field_list,overlap_bonus)



def optimize_survey(star_list,field_list,overlap_bonus):
    """
    Designs a survey that maximizes the survey score
    
    Parameter
    ---------
    star_list : list
        List of Star objects that you want to measure
    field_list : list
        List of initial fields
    overlap_bonus : float
            bonus you want to give to stars that appear in multiple fields
    Returns
    -------
    field_list : list
        Final list of Field objects
    score : float
        Score of the survey
    """
    
    N_fields = len(field_list)
    
    #Give it a list of ras, and then decs
    theta_0 = []
    
    for field in field_list:
        theta_0.append(field.ra)
        
    for field in field_list:
        theta_0.append(field.dec)
    
    field_size = field_list[0].size
        
    res = minimize(f_min,theta_0,args=([star_list,field_size,overlap_bonus]))
    
    theta = res['x']
    
    score = make_survey_and_score(theta,(star_list,field_size,overlap_bonus))
    
    ras = theta[:N_fields]
    decs = theta[N_fields:]
    
    field_list = [Field(ra,dec,field_size) for ra,dec in zip(ras,decs)]
    
    region_string_list = [f.to_region_string(star_list = star_list) for f in field_list]
    
    return field_list,region_string_list,score

def check_overlap(field_list, proposed_field):
    """
    If proposed_field overlaps with any of the fields in field_list by too much, 
    returns False
    
    Parameters
    ----------
    field_list : list
        List of Field objects
    proposed_field : Field
    
    Returns
    -------
    overlap : bool
        True if none of the Fields overlap, else False
    """
    
    field_size_deg = proposed_field.size / 60
    
    dec_seps = np.array([np.abs(f.dec - proposed_field.dec) for f in field_list])
    ra_seps = np.cos(np.pi*proposed_field.ra/180)*np.array([np.abs(f.dec - proposed_field.dec) for f in field_list])
    
    if np.any(dec_seps < field_size_deg/2) or np.any(ra_seps < field_size_deg/2):
        return False
    
    return True

def initialize_fields(N_fields,size,star_list,optimize = False, overlap = 0.1):
    """
    Makes a survey by randomly choosing fields. If any of the fields overlap by half of
    size, tries again
    
    Parameters
    ----------
    N_fields : int
        Number of fields you want
    size : float
        Size of fields in arcminutes
    star_list : list
        List of Star objects you want to observe
    optimize : bool
        If you want to run optimize_survey on the final list to tweak things a little
    overlap : float
        Overlap bonus to apply when optimizing survey
        
    Returns
    -------
    field_list : list
        List of Field objects
    
    """
    
    field_list = [Field.from_star_list(star_list=star_list,size=size)]
    
    fail_count = 0
    
    while len(field_list) < N_fields:
        
        test_field = Field.from_star_list(star_list=star_list,size=size)
        
        if check_overlap(field_list,test_field):
            
            fail_count = 0
            
            field_list.append(test_field)
            
        else:
            
            fail_count += 1
            
        if fail_count >= 100:
            
            break
    
    while len(field_list) != N_fields:
        #If we failed too many times, just add fields until we're good.
        field_list.append(Field.from_star_list(star_list=star_list,size=size))
    
    if optimize:
        field_list = optimize_survey(star_list,field_list,overlap)[0]
        
    return field_list

def tile_stars(star_list, field_size, overlap_fraction=0.05, margin_fraction=0.1):
    """
    Creates a survey by tiling a field, then deletes the fields that contain no stars
    
    Parameters
    ----------
    star_list : list
        List of `Star` objects
    field_size : float
        Size of camera in arcminutes
    overlap_fraction : float
        Fraction of the size of the camera to overlap the fields for handy image alignment
    margin_fraction : float
        Fraction of the camera to nudge over so that the first image has a margin
        
    Returns
    -------
    field_list : list
        List of `Field` objects
    """
    
    star_coords = np.array([[star.ra,star.dec] for star in star_list])
    
    #Figure out the furthest northeast point to cover
    max_ra = np.max(star_coords[:,0])
    min_ra = np.min(star_coords[:,0])
    max_dec = np.max(star_coords[:,1])
    min_dec = np.min(star_coords[:,1])
    
    #Add the margins to it
    top_edge = max_dec + (margin_fraction*field_size/60.0)
    #subtract off half the field size
    top_center = top_edge - ((field_size/60.0) / 2.0)
    
    #Same but include spherical geometry, ugh
    left_edge = max_ra + (margin_fraction*field_size/60.0)/np.cos(top_center * np.pi / 180)
    left_center = left_edge - ((field_size/60.0) / 2.0)/np.cos(top_center * np.pi / 180)
    
    ras = [left_center]
    decs = [top_center]
    
    current_ra = left_center
    current_dec = top_center
    
    while current_ra >= (min_ra - (margin_fraction*field_size/60.0)/np.cos(top_center * np.pi / 180)):
        
        current_ra -= (1.0 - overlap_fraction)*field_size/np.cos(top_center * np.pi / 180)/60.0
        ras.append(current_ra)
        
    while current_dec >= (min_dec - (margin_fraction*field_size/60.0)):
        
        current_dec -= (1.0 - overlap_fraction)*field_size/60.0
        decs.append(current_dec)
        
    temp_list = [Field(ra,dec,field_size) for ra in ras for dec in decs]
    field_list = []
    
    for field in temp_list:
        num_stars = field.count_stars_in_field(star_list=star_list)
        if num_stars != 0:
            field_list.append(field)
    
    return field_list

if __name__ == '__main__':
    
    member_hdu = fits.open('cluster_members.fits')
    member_table = Table(member_hdu[1].data)
    OB_table = member_table[member_table['SpT'] <= 20]

    OB_list = [Star(ra,dec) for ra,dec in zip(OB_table['RAJ2000'],OB_table['DEJ2000'])]
    
    N_fields = 30
    field_size = 3.0
    overlap_bonus = 0.1
    
    initial_list = initialize_fields(N_fields,field_size,OB_list)
    
    ndim, nwalkers = 2*N_fields, 10*N_fields
    p0 = []
    for i in range(nwalkers):
        
        pos_arr = []
        
        for field in initial_list:
            pos_arr.append(field.ra)
            
        for field in initial_list:
            pos_arr.append(field.dec)
        
        p0.append(pos_arr)
        
    p0 = np.array(p0)
    
    def lnprob(theta):

        ras = theta[:N_fields]
        decs = theta[N_fields:]

        field_list = [Field(ra,dec,field_size) for ra,dec in zip(ras,decs)]

        score = score_survey(OB_list,field_list,overlap_bonus)

        return np.log(score)

    sampler = mc.EnsembleSampler(nwalkers, ndim, lnprob)
    sampler.run_mcmc(p0, 1000)
    
    chain = sampler.flatchain
    
    np.savetxt('outchain.txt',chain)