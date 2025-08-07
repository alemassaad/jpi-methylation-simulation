"""Logging utilities for batch processing."""

import logging
import os
from datetime import datetime
from typing import Optional


def setup_logger(
    name: str = "batch_processor",
    log_dir: str = "logs",
    level: str = "INFO",
    console: bool = True
) -> logging.Logger:
    """
    Set up a logger with file and console handlers.
    
    Args:
        name: Logger name
        log_dir: Directory for log files
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        console: Whether to also log to console
        
    Returns:
        Configured logger instance
    """
    # Create log directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)
    
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))
    
    # Avoid duplicate handlers
    if logger.handlers:
        return logger
    
    # Create formatters
    detailed_formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    simple_formatter = logging.Formatter('%(levelname)-8s | %(message)s')
    
    # File handler with detailed format
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'batch_run_{timestamp}.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(detailed_formatter)
    logger.addHandler(file_handler)
    
    # Console handler with simple format
    if console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(simple_formatter)
        logger.addHandler(console_handler)
    
    logger.info(f"Logging initialized. Log file: {log_file}")
    
    return logger


def get_logger(name: str = "batch_processor") -> logging.Logger:
    """Get an existing logger instance."""
    return logging.getLogger(name)


class ProgressLogger:
    """Helper class for logging progress with indentation."""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.indent_level = 0
    
    def indent(self):
        """Increase indentation level."""
        self.indent_level += 1
    
    def dedent(self):
        """Decrease indentation level."""
        self.indent_level = max(0, self.indent_level - 1)
    
    def log(self, level: str, message: str):
        """Log a message with current indentation."""
        indent = "  " * self.indent_level
        prefix = "├── " if self.indent_level > 0 else ""
        self.logger.log(
            getattr(logging, level.upper()),
            f"{indent}{prefix}{message}"
        )
    
    def info(self, message: str):
        """Log info message with indentation."""
        self.log("INFO", message)
    
    def success(self, message: str, duration: Optional[float] = None):
        """Log success message with optional duration."""
        if duration:
            message = f"{message} OK ({duration:.1f}s)"
        else:
            message = f"{message} OK"
        self.log("INFO", message)
    
    def error(self, message: str):
        """Log error message with indentation."""
        self.log("ERROR", message)
    
    def warning(self, message: str):
        """Log warning message with indentation."""
        self.log("WARNING", message)