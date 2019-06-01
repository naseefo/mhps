import click
from data.defaults.params import get_default_param_values, view_default_param_values, set_default_param_values, delete_default_value


@click.group()
def default():
    pass


@default.command('view')
@click.argument('var_name',default=".", type=str)
def view_default(var_name):
    view_default_param_values(var_name.upper())

@default.command('set')
@click.argument('var_name', type=str)
@click.argument('var_values',nargs=-1)
def set_default(var_name, var_values):
    var_values_type_correction = []
    for idx, val in enumerate(var_values):
        typed_val = val
        try:
            typed_val = float(val)
        except:
            pass
        try:
            typed_val = int(val)
        except:
            pass
        var_values_type_correction.append(typed_val)
    
    if len(var_values_type_correction) == 1:
        set_default_param_values(var_name.upper(), var_values_type_correction[0])
    else:
        set_default_param_values(var_name.upper(), var_values_type_correction)

@default.command('del')
@click.argument('var_name', default=None, type=str)
def delete_default(var_name):
    if var_name is not None:
        delete_default_value(var_name.upper())