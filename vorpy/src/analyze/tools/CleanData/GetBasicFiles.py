import os
import sys
from pathlib import Path
from typing import Optional, List, Union
import logging

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))

from Data.Analyze.tools.CleanData.utils import (
    get_directory,
    safe_copy,
    ensure_directory,
    logger
)

# Configure file-specific logger
logger = logging.getLogger(__name__)

class BasicFileCollector:
    """Collects and organizes basic files from a source directory structure."""
    
    REQUIRED_FILES = [
        'balls.txt',
        'balls.pdb',
        'info.txt',
        'retaining_box.off',
        'set_atoms.pml',
        'set_balls.pml'
    ]
    
    AW_FILES = [
        'aw/aw_verts.txt',
        'aw/aw_logs.csv'
    ]
    
    POW_FILES = [
        'pow/pow_verts.txt',
        'pow/pow_logs.csv'
    ]
    
    def __init__(self, source_dir: Optional[Union[str, Path]] = None):
        """
        Initialize the collector with a source directory.
        
        Args:
            source_dir: Path to source directory. If None, will prompt for selection.
        """
        self.source_dir = Path(source_dir) if source_dir else get_directory(title="Select Source Directory")
        self.basic_data_dir = self.source_dir.parent / 'Basic_Data'
        
    def setup_basic_data_directory(self) -> None:
        """Create the Basic_Data directory and required subdirectories."""
        ensure_directory(self.basic_data_dir)
        
    def copy_files_for_subfolder(self, subfolder: Path) -> None:
        """
        Copy required files for a specific subfolder.
        
        Args:
            subfolder: Path to the subfolder to process
        """
        if subfolder.name in {'foam_data.csv', 'overlaps.csv'}:
            return
            
        target_subfolder = self.basic_data_dir / subfolder.name
        ensure_directory(target_subfolder)
        ensure_directory(target_subfolder / 'aw')
        ensure_directory(target_subfolder / 'pow')
        
        # Copy main files
        for filename in self.REQUIRED_FILES:
            safe_copy(subfolder / filename, target_subfolder / filename)
            
        # Copy AW files
        for filename in self.AW_FILES:
            if not safe_copy(subfolder / filename, target_subfolder / filename):
                # Try without subdirectory
                safe_copy(subfolder / filename.split('/')[-1], 
                         target_subfolder / filename)
                
        # Copy POW files
        for filename in self.POW_FILES:
            if not safe_copy(subfolder / filename, target_subfolder / filename):
                # Try without subdirectory
                safe_copy(subfolder / filename.split('/')[-1],
                         target_subfolder / filename)
    
    def process_all(self) -> None:
        """Process all subfolders in the source directory."""
        self.setup_basic_data_directory()
        
        subfolders = [f for f in self.source_dir.iterdir() if f.is_dir()]
        total = len(subfolders)
        
        logger.info(f"Processing {total} folders...")
        
        for i, subfolder in enumerate(subfolders, 1):
            logger.info(f"Processing folder {i}/{total}: {subfolder.name}")
            try:
                self.copy_files_for_subfolder(subfolder)
            except Exception as e:
                logger.error(f"Error processing {subfolder.name}: {str(e)}")
                continue

        logger.info("Processing complete!")

def main():
    """Main entry point for the script."""
    try:
        collector = BasicFileCollector()
        collector.process_all()
    except Exception as e:
        logger.error(f"Error during execution: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()

