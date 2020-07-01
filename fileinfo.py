'''Read paths from paths.conf, or support override in environment variables'''
from __future__ import absolute_import
PATH_CONFIG_FILE='paths.conf'

prefix='stokes-numerics-data/'

# Dict mapping path variable names to default /tmp suffix
PATH_VARS = {
    'XARPATH': prefix+'xars',
    'FRAMEPATH': prefix+'frames',
    'HKMETRICPATH': prefix+'hk',
    'CACHEPATH': prefix+'cache',
}

try:
    import configparser
except ImportError:  # python2 will raise this exception
    import ConfigParser as configparser
import os
import errno
import logging
import tempfile
cp = configparser.RawConfigParser()
logger = logging.getLogger(__name__)

# We read variables in the script directory config file first, but
# then let a config file in the current working directory override
# them (printing a warning).
config_fn_scriptdir = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),PATH_CONFIG_FILE))
config_fn_cwd = os.path.realpath(PATH_CONFIG_FILE)

cp.read([config_fn_scriptdir])
cp.read([config_fn_cwd])

if (config_fn_scriptdir != config_fn_cwd) and os.path.exists(config_fn_scriptdir) and os.path.exists(config_fn_cwd):
    logger.warning('Settings in "%s" override those in "%s" (working_dir takes priority over module_dir)' % (config_fn_cwd, config_fn_scriptdir))

for p in PATH_VARS:
    try:
        v = cp.get('paths',p)
    except configparser.Error:
        v = None
    venv = os.getenv(p)
    if venv != None:
        v = venv

    if v == None:
        v = os.path.join(tempfile.gettempdir(),PATH_VARS[p])
        logger.warning('"%s" not configured in paths.conf or environment variable; using default "%s" instead.' % (p,v))

    globals()[p] = v

    # create the directory if needed
    if not os.path.exists(v):
        logger.warning('Creating "%s"' % v)
        # avoid race condition by catching EEXIST
        # (v may have been created between check and now)
        # In python3 could use os.makedirs(...,exists_ok=True)
        try:
            os.makedirs(v)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
