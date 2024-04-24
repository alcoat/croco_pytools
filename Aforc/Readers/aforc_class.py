class aforc_class:
    def __init__(self):
        self.raw_name = {}
        self.units = {}
        self.conv_cff = {}

    def add_var(self, key, value):
        self.raw_name[key] = value

    def set_units(self, key, unit):
        self.units[key] = unit

    def set_conv(self, key, conv_cff):
        self.conv_cff[key] = conv_cff

    def get_var(self, key):
        return self.raw_name.get(key, None)

    def get_units(self, key):
        return self.units.get(key, None)

    def get_conv_cff(self, key):
        return self.conv_cff.get(key,None)

def create_class(var_info):
    variables = aforc_class()
    for i in range(len(var_info)):
        variables.add_var(var_info[i][0],var_info[i][1])
        variables.set_units(var_info[i][0],var_info[i][2])
        variables.set_conv(var_info[i][0],var_info[i][3])
    return variables



