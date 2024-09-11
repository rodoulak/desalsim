# Scaling Functions for Technology Processes
"""
    Scale-up function for technology-specific processes.
    
    This function calculates the scaled-up equipment cost or other relevant 
    parameter (C2) based on the equipment cost (C1), the original 
    mass or volume (M1), the scaled mass or volume (M2), and the specific 
    technology being used.
    
    Parameters:
    -----------
    C1 : float
        The original equipment cost or value of the parameter to be scaled.
    
    M1 : float
        The original mass or volume associated with C1.
    
    M2 : float
        The scaled mass or volume for which the new concentration or 
        parameter value (C2) will be calculated.
    
    x : float
        The scaling factor is defined based on the technology/process. 
        For these technologies, the scaling factor x is set to 0.6.
    
    Returns:
    --------
    C2 : float
        The scaled concentration or value of the parameter after applying 
        the scaling factor.
"""

def scaleup(C1,M1,M2):
    x=0.6        
    C2=C1*(M2/M1)**x
    return C2