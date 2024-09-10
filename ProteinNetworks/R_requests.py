import os
import pandas as pd
from   subprocess import Popen, PIPE


# path to directory contains all RScripts
RSCRIPTS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'RScripts')

def Check_R_packages(CRAN_packages:None, BiocManager_packages:None) -> bool:
    """
    function checks target R-packages among installed R-packages
    function makes CMD request to RScript
    :param CRAN_packages: list of CRAN packages to check
    :param BiocManager_packages: list of BiocManager packages to check
    :return:
    """
    # if CRAN_packages is None, then CRAN_packages = empty list
    CRAN_packages = CRAN_packages or []
    BiocManager_packages = BiocManager_packages or []

    s = ''
    for p in CRAN_packages:
        s += str(p) + ','  # s - temp string contained all target CRAN-package names

    if len(s) > 0:
        CRAN = '1'
        with open('temp_CRAN_packages.txt', 'w+') as f1:
            f1.write(s.rstrip(','))
            f1.write('\n')
    else:
        CRAN = '0'

    s = ''
    for p in BiocManager_packages:
        s += str(p) + ',' # s - temp string contained all target BiocManager-package names

    if len(s) > 0:
        BiocM = '1'
        with open('temp_BiocManager_packages.txt', 'w+') as f2:
            f2.write(s.rstrip(','))
            f2.write('\n')
    else:
        BiocM = '0'

    command = 'Rscript'
    path2script = os.path.join(RSCRIPTS_PATH, 'check_packages.R')

    # Variable number of args in a list
    args = [CRAN, BiocM]
    # Build subprocess command
    cmd = [command, path2script] + args
    # check_output will run the command and store to result
    p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output, error = p.communicate()

    # PRINT R CONSOLE OUTPUT (ERROR OR NOT)
    if p.returncode == 0 and len(output) > 2:
        print(f'R OUTPUT:\n {output.decode("utf8")}')
        # TEMPORARY SOLUTION
        # Check for not installed packages message:
        if output.decode("utf8").find('Packages not found') != -1:
            installing = True
        else: installing = False
    else:
        if len(error)>0: print(f'R ERROR:\n {error.decode("utf8")}')
        installing = False

    try: os.remove('temp_CRAN_packages.txt')
    except: pass
    try: os.remove('temp_BiocManager_packages.txt')
    except: pass

    return installing

def short_R_output(text):
    """
    function for shortening R output string. It takes string as parameter and return new string which
    contains only lines that starts with '[1]'
    This allows you to remove system messages, leaving only messages transmitted via 'print'
    
    :param text: text to be shortened
    :return: shortened text
    """
    df = pd.DataFrame(text.split('\n'), columns=['text'])
    new_text = df.text[df.text.apply(lambda x: True if x.find('[1]') >= 0 else False )].apply(lambda x: x.strip('[1] ').rstrip('\n')).apply(
                    lambda x: str(x) + '\n').sum()
    #isn`t working, for some reason 'contains' includs extra lines
    #new_text = df.text[df.text.str.contains('[1]')].apply(lambda x: x.strip('[1] ').rstrip('\n')).apply(
    #                lambda x: str(x) + '\n').sum()
    return new_text