(() => {
  const token = document.querySelector('meta[name="csrf-token"]')?.content || "";
  const nativeFetch = window.fetch.bind(window);
  window.fetch = (input, init = {}) => {
    const method = String(init.method || "GET").toUpperCase();
    const url = typeof input === "string" ? input : input?.url || "";
    const sameOrigin = !url || new URL(url, window.location.href).origin === window.location.origin;
    if (!sameOrigin || !["POST", "PUT", "PATCH", "DELETE"].includes(method)) {
      return nativeFetch(input, init);
    }
    const headers = new Headers(init.headers || (typeof input !== "string" ? input?.headers : undefined));
    headers.set("X-CSRF-Token", token);
    return nativeFetch(input, { ...init, headers });
  };
})();
