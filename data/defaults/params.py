
import json


def get_default_param_values():
    # Get params function
    with open("data/defaults/params.json", 'r') as f1:
            default_val = json.load(f1)
    return default_val


def view_default_param_values(var_name):
    if var_name == ".":
        default_val = get_default_param_values()
        print(default_val)
    else:
        default_val = get_default_param_values()
        try:
            print(var_name +" : " + str(default_val[var_name]))
        except:
            print("No such variable exists!")


def set_default_param_values(var_name, values):
    default_val = get_default_param_values()
    default_val.update({var_name:values})
    with open("data/defaults/params.json", 'w') as outfile:
        json.dump(default_val, outfile)
    print('Successfully added {} to default file'.format(var_name))


def delete_default_value(var_name):
    default_val = get_default_param_values()
    default_val.pop(var_name, None)
    with open("data/defaults/params.json", 'w') as outfile:
        json.dump(default_val, outfile)
    print('Successfully deleted {} from default file'.format(var_name))