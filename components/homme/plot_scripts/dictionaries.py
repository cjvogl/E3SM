# This is meant to be an all-inclusive list of possible methods
methodDict = {'KGU35-native': 5,
              'ARS232-native':  7,
              'KGS252-native':  9,
              'KGS254-native': 14,
              'KGU35': 21,
              'ARS232': 22,
              'DBM453': 23,
              'ARS222': 24,
              'ARS233': 25,
              'ARS343': 26,
              'ARS443': 27,
              'ARK324': 28,
              'ARK436': 29,
              'SSP3333b': 30,
              'SSP3333c': 31,
              'KGS232': 32,
              'KGS242': 33,
              'KGS243': 34,
              'KGS252': 35,
              'KGS254': 36}

# Line+marker styles chosen to provide certain information about the method
#   line:
#     solid      = second order
#     dashed     = third order
#     dot-dashed = fourth order
#   marker:
#     circle   = 0 implicit stages
#     x        = 2 implicit stages
#     triangle = 3 implicit stages
#     square   = 4 implicit stages
#     pentagon = 5 implicit stages
lineStyleDict = {'KGU35-native': '--o',
                 'ARS232-native': '-x',
                 'KGS252-native': '-x',
                 'KGS254-native': '-s',
                 'KGU35': '--o',
                 'ARS222': '-x',
                 'ARS232': '-x',
                 'KGS232': '-x',
                 'KGS242': '-x',
                 'KGS243': '-^',
                 'KGS252': '-x',
                 'KGS254': '-s',
                 'KGU35': '--o',
                 'SSP3333b': '--x',
                 'SSP3333c': '--x',
                 'ARS443': '--s',
                 'DBM453': '--s',
                 'ARS233': '--x',
                 'ARS343': '--^',
                 'ARK324': '--^',
                 'ARK436': '-.p'}

colorDict = {'KGU35-native': 'tab:blue',
             'KGU35': 'tab:orange',
             'ARS232-native': 'tab:blue',
             'KGS252-native': 'tab:orange',
             'ARS232': 'tab:green',
             'ARS222': 'tab:red',
             'KGS232': 'tab:purple',
             'KGS242': 'tab:brown',
             'KGS252': 'tab:pink',
             'ARS233': 'tab:gray',
             'SSP3333b': 'tab:olive',
             'SSP3333c': 'tab:cyan',
             'KGS243': 'tab:green',
             'ARS343': 'tab:red',
             'ARK324': 'tab:purple',
             'KGS254-native': 'tab:brown',
             'KGS254': 'tab:pink',
             'ARS443': 'tab:gray',
             'DBM453': 'tab:olive',
             'ARK436': 'tab:cyan'}
