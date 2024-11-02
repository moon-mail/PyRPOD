import ast
import sys

def get_imports(filepath):
    """Extract all imported modules, functions, and aliases from a Python file."""
    with open(filepath, "r") as file:
        tree = ast.parse(file.read())
    
    imports = {}
    
    for node in ast.walk(tree):
        # For regular imports (e.g., `import module as alias`)
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports[alias.asname or alias.name.split('.')[0]] = alias.name.split('.')[0]
                
        # For from-imports (e.g., `from module import func as alias`)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                for alias in node.names:
                    full_name = f"{node.module}.{alias.name}"
                    imports[alias.asname or alias.name] = node.module.split('.')[0]
    
    return imports

def find_used_imports(filepath, imports):
    """Check if imports are used in the Python file."""
    used_imports = set()

    with open(filepath, "r") as file:
        lines = file.readlines()
    
    for line in lines:
        # Ignore lines that are just import statements
        if line.strip().startswith("import") or line.strip().startswith("from"):
            continue

        # Check if any imported module, function, or alias is used in each line
        for alias, module in imports.items():
            if alias in line:
                used_imports.add(module)
    
    return used_imports

def find_unused_imports(filepath):
    """Identify unused imports in a Python file."""
    imports = get_imports(filepath)
    used_imports = find_used_imports(filepath, imports)

    # Find unused imports by comparing all imports with used imports
    unused_imports = set(imports.values()) - used_imports

    # Display results
    if unused_imports:
        print("Unused imports found:")
        for imp in unused_imports:
            print(f"- {imp}")
    else:
        print("No unused imports found. All imports are used.")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_unused_imports.py <path_to_python_file>")
    else:
        filepath = sys.argv[1]
        find_unused_imports(filepath)
