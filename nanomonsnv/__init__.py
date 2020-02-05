#! /usr/bin/env python

from .parser import create_parser

def main():
    parser = create_parser()
    if not hasattr(parser.parse_args(), 'func'):
        parser.error('too few arguments')
    args = parser.parse_args()
    args.func(args)
    
    