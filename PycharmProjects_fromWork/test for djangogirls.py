import requests

url = 'http://127.0.0.1:8000/hello/'
headers = {'Authorization': 'Token ae0966b1439e8d731140c8b6278da275fcd56cf0'}
r = requests.get(url, headers=headers)
print(r.content)
url_2 = 'http://127.0.0.1:8000/api-token-auth/'
user_data = {'username':'maloseva', 'password':'astra030691'}
r_post = requests.post(url, data=user_data)
print(r_post)