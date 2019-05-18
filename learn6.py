
import logging
import boto3
from botocore.exceptions import ClientError


def create_bucket(bucket_name):
    """ Create an Amazon S3 bucket

    :param bucket_name: Unique string name
    :return: True if bucket is created, else False
    """

    s3 = boto3.client('s3', aws_access_key_id="AKIA3N2CFDTP3XCRNHUB", aws_secret_access_key= "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", region_name='us-east-2b')
    try:
        s3.create_bucket(Bucket=bucket_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True

create_bucket('mhps-iitd-asd3')