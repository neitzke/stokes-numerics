'''Get the head revision hash if we're in a git repo and GitPython is installed'''
import os
import logging
logger = logging.getLogger(__name__)
try:
    import git
    try:
        repo = git.Repo(os.path.dirname(os.path.realpath(__file__)),
                        search_parent_directories=True)
        codeversion = repo.head.object.hexsha
        if repo.is_dirty():
            codeversion += "-dirty"
    except git.InvalidGitRepositoryError:
        logger.warning('Not running in a git repository; output will not be tagged with commit hash')
        # This is not a git repository
        codeversion = "norepo"
except ImportError:
    logger.warning('GitPython is not available; output will not be tagged with a commit hash')
    # We don't know if it is a git repository
    codeversion = "unknown"
