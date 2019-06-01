
import logging
import boto3
from botocore.exceptions import ClientError


def create_bucket(bucket_name):
    """ Create an Amazon S3 bucket

    :param bucket_name: Unique string name
    :return: True if bucket is created, else False
    """

    
    try:
        s3.create_bucket(Bucket=bucket_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True

create_bucket('mhps-iitd-asd3')