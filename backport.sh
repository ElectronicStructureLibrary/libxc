#!/usr/bin/env bash

# From BerkeleyGW r5128, references to svninfo.f90 removed.

# BerkeleyGW, Copyright (c) 2011, The Regents of the University of
# California, through Lawrence Berkeley National Laboratory (subject to
# receipt of any required approvals from the U.S. Dept. of Energy).
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# (1) Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 
# (2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 
# (3) Neither the name of the University of California, Lawrence
# Berkeley National Laboratory, U.S. Dept. of Energy nor the names of
# its contributors may be used to endorse or promote products derived
# from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# You are under no obligation whatsoever to provide any bug fixes,
# patches, or upgrades to the features, functionality or performance of
# the source code ("Enhancements") to anyone; however, if you choose to
# make your Enhancements available either publicly, or directly to
# Lawrence Berkeley National Laboratory, without imposing a separate
# written license agreement for such Enhancements, then you hereby grant
# the following license: a  non-exclusive, royalty-free perpetual
# license to install, use, modify, prepare derivative works, incorporate
# into other computer software, distribute, and sublicense such
# enhancements or derivative works thereof, in binary and source code
# form.

# David Strubbe, March 2012, UC Berkeley
# Script to automagically backport a commit from trunk to a branch.
# Usage (from checked-out version of branch): sh backport.sh revision#
# Has been tested for svn versions 1.6, 1.7, 1.8

if [ $# -ne 1 ]; then
    echo "Usage: $0 revision"
    exit
fi

# Without this, svn 1.8 may give an error such as:
#svn: E195020: Cannot merge into mixed-revision working copy [5751:5752]; try updating first
svn update

#Sample output from svn info with subversion 1.6
#Path: .
#URL: http://www.tddft.org/svn/octopus/trunk
#Repository Root: http://www.tddft.org/svn/octopus
#Repository UUID: 53d8e6fd-f412-0410-a60b-8e56063eb7bf
#Revision: 10120
#Node Kind: directory
#Schedule: normal
#Last Changed Author: xavier
#Last Changed Rev: 10120
#Last Changed Date: 2013-02-25 09:36:09 -0500 (Mon, 25 Feb 2013)

#Sample output from svn info with subversion 1.7 -- as above, except after Path comes also:
#Working Copy Root Path: /Users/dstrubbe/Software/octopus

#Sample output from svn info with subversion 1.8 -- as above, except after URL comes also: 
#Relative URL: ^/trunk 

# e.g. https://civet.berkeley.edu/svn/BerkeleyGW
Root=`svn info | grep "Repository Root" | awk '{print $3}'`
# e.g. https://civet.berkeley.edu/svn/BerkeleyGW/branches/1.0.x
URL=`svn info | grep URL | awk '{print $2; exit}'`
# e.g. branches/1.0.x
Branch=${URL#$Root/}

echo "svn merge -c $1 $Root/trunk ."
svn merge -c $1 $Root/trunk . --accept postpone --quiet

# We will accumulate annoying and useless merge tracking properties otherwise
svn revert . --quiet

echo
echo "===== svn status ====="
svn status

echo "Backport of $1 to $Branch:" > svn-commit.tmp
svn log $Root -r $1 >> svn-commit.tmp

echo
echo "===== log message ====="
cat svn-commit.tmp

echo
echo "Ok to commit? [y/n]"
read response

if [ "$response" = "y" ]; then
# remove temporary file, but only if commit succeeds
  svn commit --file=svn-commit.tmp && rm -f svn-commit.tmp
else
  echo "Not committing."
fi
