#!/bin/bash

SOURCEDIR=mmayers@garibaldi:"/gpfs/home/mmayers/proteomics/N15_mouse_studies/"

rsync -av --exclude="*UC13" --include="*/" --include="DTASelect-filter*.txt" --exclude="*" --prune-empty-dirs $SOURCEDIR .

find . -type d -exec chmod 755 {} +
