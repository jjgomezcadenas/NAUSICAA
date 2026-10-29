"""Run the viewer integration checks against the current direct-collected inputs."""
import sys
from check_direct_inputs import main

if __name__=='__main__':
    # Preserve the earlier positional URL command.
    if len(sys.argv)>1 and not sys.argv[1].startswith('--'):
        sys.argv.insert(1,'--url')
    main()
