import sys
from src.cli.poscarkit import main as main_cli
from src.cli.poscarkit_interact import main as main_interact

if __name__ == "__main__":
    # Check if there are any arguments
    if len(sys.argv) > 1:
        # If there are arguments, use the CLI mode
        main_cli()
    else:
        # If there are no arguments, use the interactive mode
        main_interact()

