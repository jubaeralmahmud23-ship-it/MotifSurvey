import argparse,json,sys
from pathlib import Path
from .core import estimate

def main():
    parser=argparse.ArgumentParser(description='Known independent-node design: motif totals and covariance')
    parser.add_argument('input',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:
        if args.output and args.output.exists():raise ValueError('Output already exists')
        payload=json.loads(args.input.read_text(encoding='utf-8-sig'))
        result=estimate(**payload)
        text=json.dumps(result,indent=2,allow_nan=False)+'\n'
        if args.output:
            with args.output.open('x',encoding='utf-8') as f:f.write(text)
        else:sys.stdout.write(text)
    except (OSError,ValueError,TypeError) as exc:
        parser.exit(2,f'motifsurvey: {exc}\n')

if __name__=='__main__':main()
