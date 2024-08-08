import argparse
from model import process

def main():
    print('parsing...')
    parser = argparse.ArgumentParser(description='register and process your histological images.')
    parser.add_argument('--dir', type=str, required=True, help='Base directory')
    parser.add_argument('--structures', type=str, required=True, help='Comma-separated list of structure acronyms')
    args = parser.parse_args()

    list = [acronym.strip() for acronym in args.structures.split(',')]

    process(args.dir, list)
    
if __name__ == "__main__":
    main()