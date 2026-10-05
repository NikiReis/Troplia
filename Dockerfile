FROM python:3.12-slim

WORKDIR /app

COPY QuantyPep/requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY QuantyPep/app.py .
COPY QuantyPep/troplia_helipep.py .
COPY QuantyPep/static/ ./static/

CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8080}"]
