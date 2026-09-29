/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  // Expose Railway API URL to server-side API routes only (not browser)
  serverRuntimeConfig: {
    railwayApiUrl: process.env.RAILWAY_API_URL || "http://localhost:8000",
  },
};

export default nextConfig;
