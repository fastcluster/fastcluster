project('fastcluster', 'cpp', version: '1.3.0')

py = import('python').find_installation(pure: false)

_python_script = '''
import os, numpy
d = numpy.get_include()
try:
  d = os.path.relpath(d)
except:
  pass
print(d)
'''

_inc_numpy = run_command(py, ['-c', _python_script], check: true).stdout().strip()

py.install_sources('fastcluster.py')

py.extension_module(
    '_fastcluster',
    'src/fastcluster_python.cpp',
    include_directories: [_inc_numpy],
    install: true,
)


