let errorHandler = null;

export function configureApi({ onError } = {}) {
  errorHandler = typeof onError === "function" ? onError : null;
}

export async function requestJson(url, options = {}) {
  const { silentError = false, ...fetchOptions } = options;
  const response = await fetch(url, fetchOptions);
  const data = await response.json();
  if (!response.ok) {
    if (!silentError) errorHandler?.(data.error || "请求失败", "error");
    throw new Error(data.error || "请求失败");
  }
  return data;
}
