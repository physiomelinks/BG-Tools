#===============================================================================
#
#  CellDL and bondgraph tools
#
#  Copyright (c) 2020 - 2025 David Brooks
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#===============================================================================

import lxml.etree as etree

#===============================================================================
#===============================================================================

class Issue(Exception):
    def __init__(self, reason: str):
        super().__init__(reason)
        self.__reason = reason

    @property
    def reason(self):
        return self.__reason

def make_issue(e: Exception) -> Issue:
#=====================================
    if isinstance(e, Issue):
        return e
    issue = Issue(str(e))
    issue.__traceback__ = e.__traceback__
    return issue

#===============================================================================
#===============================================================================

"""
Generate URIs for lxml.etree.
"""
class XMLNamespace:
    def __init__(self, ns: str):
        self.__ns = ns

    def __str__(self):
        return self.__ns

    def __call__(self, attr: str='') -> str:
        return f'{{{self.__ns}}}{attr}'

    def __getattr__(self, attr: str) -> str:
        return f'{{{self.__ns}}}{attr}'

    def url(self, suffix: str) -> str:
        return f'{self.__ns}{suffix}'

#===============================================================================

def etree_from_string(xml: str) -> etree.Element:
    parser = etree.XMLParser(remove_blank_text=True)
    try:
        return etree.fromstring(xml, parser)
    except etree.XMLSyntaxError as error:
        raise ValueError(error)

#===============================================================================

def string_to_list(string: str) -> list[int]:
#============================================
    return list(string.encode('utf-8'))

#===============================================================================
#===============================================================================
