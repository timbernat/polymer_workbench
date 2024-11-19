'''Dynamically compile and store all imports occurring within polymerist'''

import sys
import pandas as pd

import polymerist as ps
from polymerist.genutils.importutils.pyimports import extract_imports_from_module


infos = extract_imports_from_module(ps)
df = pd.DataFrame.from_records([info.__dict__ for info in infos])
df.to_csv('imports.csv')

# separating relative from non-relative imports
nonrel = [info for info in infos if not info.is_relative and info.parent_module is None]
print(f'{len(nonrel)} non-relative imports found')

imported_names = set(info.object_name for info in nonrel)
imported_names

registered_builtins = set(sys.builtin_module_names)
registered_stdlibs = set(sys.stdlib_module_names)

nb_imports = imported_names - registered_builtins - registered_stdlibs
print(f'Unique non-relative imports found are:\n {nb_imports}')