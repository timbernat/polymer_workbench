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
            working_line = 0 if mdn is None else mdn.end_lineno

            i = 0
            found_nonempty = False
            for line in lines[working_line:]:
                if line != '\n':
                    found_nonempty = True
                    break
                i += 1
            print(path, working_line + i - num_lines, found_nonempty)

            if i < 1 and found_nonempty:
                lines.insert(working_line, '\n') # if no gap is present, insert a newline
                i += 1
            if i > 1:
                for _ in range(i - 1):
                    lines.pop(working_line)
                    i -= 1
            
            # 2) Make necessary insertions into text
            lineidx_nonempty = min(working_line + i, num_lines)
            insertions = [] # will be in reverse order, to avoid needing to keep track of highest insertion index

            email_nodes = [node for node in extract_dunder_tags(syntax_tree) if node.targets[0].id == '__email__']
            email_node = email_nodes[0] if email_nodes else None
            if email_node is None:
                insertions.append("__email__ = 'timotej.bernat@colorado.edu'\n")

            author_nodes = [node for node in extract_dunder_tags(syntax_tree) if node.targets[0].id == '__author__']
            author_node = author_nodes[0] if author_nodes else None
            if author_node is None:
                insertions.append("__author__ = 'Timotej Bernat'\n")
            else:
                working_line = author_node.lineno

            if insertions:
                insertions.append('\n')
            if not found_nonempty:
                insertions.append('\n')

            for added_line in insertions:
                lines.insert(working_line, added_line)

            # 3) Rewrite file
            with open(path, 'w') as file:
                lines = file.write(''.join(lines))

if __name__ == '__main__':
    main()