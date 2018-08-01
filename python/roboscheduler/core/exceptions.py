# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date:   2017-12-05 12:01:21
# @Last modified by:   Brian Cherinka
# @Last Modified time: 2017-12-05 12:19:32

from __future__ import print_function, division, absolute_import


class RoboschedulerError(Exception):
    """A custom core Roboscheduler exception"""

    def __init__(self, message=None):

        message = 'There has been an error' \
            if not message else message

        super(RoboschedulerError, self).__init__(message)


class RoboschedulerNotImplemented(RoboschedulerError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):

        message = 'This feature is not implemented yet.' \
            if not message else message

        super(RoboschedulerNotImplemented, self).__init__(message)


class RoboschedulerAPIError(RoboschedulerError):
    """A custom exception for API errors"""

    def __init__(self, message=None):
        if not message:
            message = 'Error with Http Response from Roboscheduler API'
        else:
            message = 'Http response error from Roboscheduler API. {0}'.format(message)

        super(RoboschedulerAPIError, self).__init__(message)


class RoboschedulerApiAuthError(RoboschedulerAPIError):
    """A custom exception for API authentication errors"""
    pass


class RoboschedulerMissingDependency(RoboschedulerError):
    """A custom exception for missing dependencies."""
    pass


class RoboschedulerWarning(Warning):
    """Base warning for Roboscheduler."""


class RoboschedulerUserWarning(UserWarning, RoboschedulerWarning):
    """The primary warning class."""
    pass


class RoboschedulerSkippedTestWarning(RoboschedulerUserWarning):
    """A warning for when a test is skipped."""
    pass


class RoboschedulerDeprecationWarning(RoboschedulerUserWarning):
    """A warning for deprecated features."""
    pass
