from acheron.download import are_equal_mic

def to_float(MIC):
    """
    Ensures input can be compared to float, in case input is
    string with inequalities
    """
    MIC=str(MIC)
    for equality_smb in ['=','<','>']:
        MIC = MIC.replace(equality_smb,'')
    return float(MIC)

def is_equiv(a,b):
    """
    Determins if 2 MIC's are functionally equivalent
    Allows for 1% error due to float casting

    See also acheron.download.are_equal_mic
    """
    allowed_error = 0.1
    if abs(1-(a/b)) < allowed_error:
        return True
    else:
        return False


def to_resistance_type(value, bp):
    """
    Compares value to breakpoints, returns S I or R
    Allows for 1% error do to float casting
    """
    susc_bp = bp[0]
    int_bp = bp[1]
    res_bp = bp[2]

    if is_equiv(value, susc_bp) or value < susc_bp:
        return 'S'
    elif any([is_equiv(value, i) for i in int_bp]):
        return 'I'
    elif is_equiv(value, res_bp) or value > res_bp:
        return 'R'
    else:
        raise Exception("MIC value {} not found in breakpoints {}".format(value,bp))

def find_error_type(predicted, actual, abx):
    """
    The following AST breakpoints are based on:

    CLSI M100-ED30:
    2020 Performance Standards for Antimicrobial Susceptibility Testing,
    30th Edition

    Drugs are listed in a <=S / ==I / >=R format
    May include multiple, or zero, intermediate concentrations

    NOTE:
    This is for salmonella, and will change for
    QUINOLONES AND FLUOROQUINOLONES

    Not available: Streptomycin, nalidixic acid, ceftiofur
    """

    breakpoints = {}

    """
    Penicillins
    """
    # Ampicillin
    breakpoints['AMP'] = [8,[16],32]

    """
    Beta-lactam combo's
    """
    # Amoxicillin-clavulanate
    breakpoints['AMC'] = [8,[16],32]

    """
    Cephems (cephalosporins I -IV)
    """
    # Ceftriaxone
    breakpoints['CRO'] = [1,[2],4]
    # Cefoxitin
    breakpoints['FOX'] = [8,[16],32]

    """
    Aminoglycosides
    """
    # Gentamicin
    breakpoints['GEN'] = [4,[8],16]
    # Kanamycin
    breakpoints['KAN'] = [16,[32],64]
    # Streptomycin CLSI DOES NOT LIST MIC BREAKPOINTS FOR THIS ABX
    #breakpoints['STR'] = [,[],]

    """
    Macrolides
    """
    # Azithromycin
    breakpoints['AZM'] = [16,[],32]

    """
    Tetracyclines
    """
    # Tetracycline
    breakpoints['TET'] = [4,[8],16]

    """
    QUINOLONES AND FLUOROQUINOLONES for salmonella
    """
    # Ciprofloxacin
    breakpoints['CIP'] = [0.06,[0.12,0.25,0.5],1]

    """
    FOLATE PATHWAY ANTAGONISTS
    """
    # Trimethoprim-sulfamethoxazole
    breakpoints['SXT'] = [2,[],4]
    # sulfonamides (using sulfisoxazole)
    breakpoints['FIS'] = [256,[],512]

    """
    PHENICOLS
    """
    # Chloramphenicol
    breakpoints['CHL'] = [8,[16],32]

    pred = to_float(predicted)
    act = to_float(actual)

    pred_class = to_resistance_type(pred, breakpoints[abx])
    act_class = to_resistance_type(act, breakpoints[abx])

    if is_equiv(pred, act):
        return "Correct"
    elif pred_class == 'R' and act_class == 'S':
        return "Major Error"
    elif pred_class == 'S' and act_class == 'R':
        return "Very Major Error"
    else:
        return "Non Major Error"
