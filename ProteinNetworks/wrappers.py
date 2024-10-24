import functools
import pandas as pd
import os

def display_df(df):
    """
    Function for displaying DataFrames (df). If IPNB is used, df will display with common IPNB function 'Display', else:
    it will display by print() function

    Parameters
    ----------
    df : DataFrame
        DataFrame to display

    Returns
    -------
    None
    """
    try: display(df)
    except: print(df)

def Check_Value(val:[str, float, int], valid_values:set, valname:str, message='Wrong value123'):
    """
    Function for check correctness of input value

    Parameters
    ----------
    val : str, float, int
        input value
    valid_values : set
        set of valid values
    valname : str
        group name of valid_values set. Or name of val variable
    message : str
        Error message

    Returns
    -------
    None
    """
    if val not in valid_values:
        if message == 'Wrong value123':
            message = f'Wrong value of "{valname}" variable! Choose one of {valid_values}'
        raise Exception(message)

def Check_type(val, valid_types:tuple, valname:str, message='Wrong value123'):
    
    """
    Function for check correctness of input type

    Parameters
    ----------
    val : any
        input value
    valid_types : tuple
        tuple of valid types
    valname : str
        group name of valid_values set. Or name of val variable
    message : str
        optional. Error message

    Returns
    -------
    None
    """
    if not isinstance(val, valid_types):
        if message == 'Wrong value123':
            message = f'Wrong type of "{valname}" variable! Valid types are one of {valid_types}'
        raise Exception(message)
    
    
def Check_kwargs(kwargs, valid_params:set):
    """
    Function for check correctness of input kwargs
    Parameters
    ----------
    kwargs : dict
        dict of kwargs
    Returns
    -------
    None
    """
    for param in kwargs.keys():
        if param not in valid_params:
            raise ValueError(f' {param} parameter not valid. List of valid parameters: {valid_params}')

def print_downline(line_length:int=40):
    """
    Function prints line for titler decorator

    Parameters
    ----------
    line_length : int
        length of line (number of '_' symbols)

    Returns
    -------
        None
    """
    line = '_'*line_length
    print(f'{line}\n\n')

def print_upline(title:str, line_length:int=40):
    """
    function prints line and adds title for titler decorator

    Parameters
    ----------
    line_length : int
        length of line (number of '_' symbols)

    Returns
    -------
        None
    """
    line = '_'*line_length
    print(f'\t{title}\n{line}')

def titler(title: str, line_length=40):
    """
    Decorator added Title and edges of Paragraph

    Parameters
    ----------
    title : str
        title
    line_length : int
        length of line (number of '_' symbols)

    Returns
    -------
    callable
        decorated function
    """
    def titler_decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            print_upline(title, line_length)
            original_result = func(*args, **kwargs)
            print_downline(line_length)

            return original_result
        return wrapper
    return titler_decorator

def save_table(table, name, saveformat='xlsx', index:bool = True):
    """
    Function for saving DataFrame tables
    
    Parameters
    ----------
    table : DataFrame
        DataFrame table to save
    name : str
        name of file
    saveformat : str
        format of saving file: 'xlsx' or 'csv'
    index : bool
        show indexes in saved table?
    
    Returns
    -------
    None
    """

    Check_Value(saveformat, {'csv', 'xlsx'}, 'saveformat')
    Check_type(table, pd.DataFrame, 'table')
    try:
        if saveformat == 'xlsx':
            if name[-5:] != '.xlsx' and name[-4:] != '.xls':
                name += '.xlsx'
            table.to_excel(name, index=index)
        elif saveformat == 'csv':
            if name[-4:] != '.csv':
                name += '.csv'
            table.to_csv(name, index=index, header=True)
        print(f'File {name} successfully saved in {os.path.abspath(os.getcwd())}\n')
    except PermissionError:
        print('Permission Denied Error: Access is denied. Close file if it`s open and try again')
    except:
        print('Saving file isn`t complete. If you rewrite file, close it and try again')
        
        
def create_subframe_by_names(df, column: str, names: [list, tuple, set], add: str = 'first'):
    """
    Function finds rows in original dataset and returns sub-dataframe including input names in selected column

    Parameters
    ----------
    df : pd.DataFrame
        target DataFrame
    column : str
        the selected column in which names will be searched
    names : list of str
        list of target names whose records need to be found in the table
    add : str
        ['first', 'last', 'all'] parameter of adding found rows.
        'first' - add only the first entry
        'last' - add only the last entry
        'all' - add all entries

    Returns
    -------
    pd.DataFrame
        sub-dataframe including input names in selected column
    """
    
    Check_Value(add, {'first', 'last', 'all'}, add)

    def add_all(table, rows):
        return pd.concat([table, rows])

    def add_first(table, rows):
        table.loc[len(table)] = rows.iloc[0]
        return table

    def add_last(table, rows):
        table.loc[len(table)] = rows.iloc[-1]
        return table

    adding_method = {'first': add_first,
                        'last': add_last,
                        'all': add_all}

    new_df = pd.DataFrame(columns=df.columns)
    not_found_names = []
    for name in names:
        rows = df[df[column] == name]
        if len(rows) > 0:
            new_df = adding_method[add](new_df, rows)
        else: not_found_names.append(name)
    print(f'{len(not_found_names)} names were not found in the dataframe:\n')
    print('[', end='')
    print(*not_found_names, sep=', ', end='')
    print(']')

    return new_df
