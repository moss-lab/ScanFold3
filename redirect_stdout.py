import sys, os
import logging
from instant import inline

print 'This is printed from python to stdout'
stdout = os.fdopen(os.dup(sys.stdout.fileno()), 'w')
stderr = os.fdopen(os.dup(sys.stderr.fileno()), 'w')

logging.basicConfig(stream=stderr, level=logging.DEBUG)

redirect = inline("""                                                                                                                    
void redirect(void) {                                                                                                                    
    freopen("my_stdout.txt", "w", stdout);                                                                                               
    freopen("my_stderr.txt", "w", stderr);                                                                                               
}                                                                                                                                        
""")
redirect()

cout = inline("""                                                                                                                        
void cout(void) {                                                                                                                        
    std::cout << "This is written from C++ to my_stdout.txt" << std::endl;                                                               
    std::cerr << "This is written from C++ to my_stderr.txt" << std::endl;                                                               
}                                                                                                                                        
""")
cout()

print 'This is written from python to my_stdout.txt'

stdout.write('This is printed from python to stdout\n')
stderr.write('This is printed from python to stderr\n')
logging.info('This is printed to stderr from python using logging')

"""
expected output:
$ python test.py
This is printed from python to stdout
This is printed from python to stdout
This is printed from python to stderr
INFO:root:This is printed to stderr from python using logging
$ cat my_stdout.txt 
This is written from C++ to my_stdout.txt
This is written from python to my_stdout.txt
$ cat my_stderr.txt 
This is written from C++ to my_stderr.txt
"""
