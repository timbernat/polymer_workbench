'''Script for authomatically inserting __author__ and __email__ info into all polymerist module files'''

from typing import Optional, Union

import ast
from pathlib import Path


def extract_module_docstring(module_node : ast.Module) -> Optional[ast.Expr]:
    '''
    Locates and returns a module-level docstring Expr node
    Returns None if one is not found_nonempty
    '''
    if len(module_node.body) == 0:
        return

    module_docstring_node = module_node.body[0] # expect this to come first
    if not isinstance(module_docstring_node, ast.Expr) or module_docstring_node.value.value.strip('\n') != ast.get_docstring(module_node):
        module_docstring_node = None

    return module_docstring_node

def extract_first_import(module_node : ast.Module) -> Union[ast.Import, ast.ImportFrom]:
    '''
    Locates and returns the first import statement node in a file
    Return None if no imports are present
    '''
    import_node = None
    for syntax_node in module_node.body:
        if isinstance(syntax_node, (ast.Import, ast.ImportFrom)):
            import_node = syntax_node
            break
    
    return import_node

def extract_dunder_tags(module_node : ast.Module) -> list[ast.Assign]:
    '''
    Locates and returns list of all double-underscore ("under") module attributes
    Returns empty list if none are present
    '''
    return [
        syntax_node
            for syntax_node in module_node.body
                if isinstance(syntax_node, ast.Assign) and syntax_node.targets[0].id.startswith('__')
    ]

def main():
    for path in Path('polymerist/polymerist').glob('**/*.py'):
        with open(path, 'r') as file:
            # 0) Load module text
            with open(path, 'r') as file:
                lines = file.readlines()
                num_lines = len(lines)
                syntax_tree = ast.parse(''.join(lines))
        
            # 1) Balance newlines between initial line and next non-empty line
            mdn = extract_module_docstring(syntax_tree)
            working_line_idx = -1 if mdn is None else mdn.end_lineno - 1
            if working_line_idx != -1:
                docstring_line = lines[working_line_idx]
                if not docstring_line.endswith('\n'):
                    lines[working_line_idx] += '\n'
            
            # 2) Make necessary insertions into text
            insertions = []

            author_nodes = [node for node in extract_dunder_tags(syntax_tree) if node.targets[0].id == '__author__']
            author_node = author_nodes[0] if author_nodes else None
            if author_node is None:
                insertions.append("__author__ = 'Timotej Bernat'\n")
            else:
                working_line_idx = author_node.lineno # 1-index accounts for offset here

            email_nodes = [node for node in extract_dunder_tags(syntax_tree) if node.targets[0].id == '__email__']
            email_node = email_nodes[0] if email_nodes else None
            if email_node is None:
                insertions.append("__email__ = 'timotej.bernat@colorado.edu'\n")
                
            if lines[working_line_idx] != '\n':
                insertions.append('\n')

            insertions = insertions[::-1]  # reverse order to avoid needing to keep track of highest insertion index
            if insertions:
                insertions.append('\n') # prepend newline if insertions are to be made (this will be written last)

            for added_line in insertions:
                lines.insert(working_line_idx + 1, added_line)

            print(path)
            print(insertions)
            print(lines[:10])

            # 3) Rewrite file
            with open(path, 'w') as file:
                lines = file.write(''.join(lines))

if __name__ == '__main__':
    main()