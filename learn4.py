

import os
import sys
import boto3

local_directory = 'results'
destination = 'System-a-001'

s3 = boto3.client('s3', aws_access_key_id='AKIA3N2CFDTP6OZHB2NF', aws_secret_access_key= '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')


bucket = s3.Bucket("phdresults")
bucket_prefix="System-a-001"
objs = bucket.objects.filter(Prefix =bucket_prefix)