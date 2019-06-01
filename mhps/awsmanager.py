
import os
import sys
import boto3
import logging
from botocore.exceptions import ClientError
import pandas as pd



def create_bucket(bucket_name, access_id, access_secret):
 
    print('I am in create bucket')
    data = pd.read_csv("accesskeys.csv") 
    s3 = boto3.client('s3', aws_access_key_id=access_id, aws_secret_access_key= access_secret, region_name='us-east-1')
    try:
        s3.create_bucket(Bucket=bucket_name)
        print('Bucket creation successful!')
        return s3
    except ClientError as e:
        print('Bucket creation failed!')
        logging.error(e)
        exit()

def upload(foldername, access_id, access_secret):
    # print('I am here')
    client = create_bucket(foldername, access_id, access_secret)
    print('I have reached back to upload')
    # get an access token, local (from) directory, and S3 (to) directory
    # from the command-line
    local_directory = os.path.join('results', foldername)
    bucket = foldername
    destination=''

    # client = boto3.client('s3', aws_access_key_id=access_id, aws_secret_access_key= access_secret, region_name='us-east-1')

    # enumerate local files recursively
    for root, dirs, files in os.walk(local_directory):

        for filename in files:

            # construct the full local path
            local_path = os.path.join(root, filename)
          

            # construct the full Dropbox path
            relative_path = os.path.relpath(local_path, local_directory)
         
            s3_path = os.path.join(destination, relative_path)
         

            # relative_path = os.path.relpath(os.path.join(root, filename))

            print('Searching "%s" in "%s"' % (s3_path, bucket))
            try:
                client.head_object(Bucket=bucket, Key=s3_path)
                print("Path found on S3! Skipping %s..." % s3_path)

                # try:
                    # client.delete_object(Bucket=bucket, Key=s3_path)
                # except:
                    # print "Unable to delete %s..." % s3_path
            except:

                print("Uploading %s..." % s3_path)
                client.upload_file(local_path, bucket, s3_path)
                print('Upload successful!')